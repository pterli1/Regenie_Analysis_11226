"""DNAnexus job and file helpers, wrapping dxpy.

Calls to the DNAnexus API are made one at a time, except job termination.
"""

from __future__ import annotations

import asyncio
import concurrent.futures
import math
import re
from collections import Counter
from pathlib import Path

import dxpy
import dxpy.api
from loguru import logger

from fedcat.core.logging import NCPU, Permits

_DXPY_CONFIGURED: dict[tuple[str, str], bool] = {}
_DXAPI_SEMAPHORE = asyncio.BoundedSemaphore(1)
_DXAPI_TERMINATE_SEMAPHORE = asyncio.Semaphore(8)
_DESCRIBE_BATCH_SIZE = 500
_DXPY_SETUP_LOCK = asyncio.Lock()

UKB_APPLET_PROJECT = "project-GbZ9g7QJVJzZVQj45j6xgxZ8"
"""Default DNAnexus project holding the applets."""


async def setup_dxpy(project_requirements: dict[str, str]) -> None:
    """Check the dxpy login and the access level to each project, once per project and level.

    Args:
        project_requirements: Minimum permission level per project ID: `"VIEW"`,
            `"UPLOAD"`, `"CONTRIBUTE"` or `"ADMINISTER"`.

    Raises:
        PermissionError: If access to a project is below its level.
        ValueError: If a level is not valid.
        dxpy.exceptions.DXError: If not logged in, or a project cannot be described.
    """

    if all(_DXPY_CONFIGURED.get((pid, level), False) for pid, level in project_requirements.items()):
        logger.trace("dxpy already configured, skipping setup.")
        return

    async with _DXPY_SETUP_LOCK:
        if all(_DXPY_CONFIGURED.get((pid, level), False) for pid, level in project_requirements.items()):
            logger.trace("dxpy already configured, skipping setup.")
            return

        permission_levels = {"NONE": 0, "VIEW": 1, "UPLOAD": 2, "CONTRIBUTE": 3, "ADMINISTER": 4}

        try:
            user = await asyncio.to_thread(dxpy.whoami)
            logger.info(f"Authenticated with DNAnexus as {user}")
        except Exception:
            logger.error("dxpy not set up, make sure to login on command line")
            raise

        for project_id, min_level in project_requirements.items():
            logger.info(f"Checking project {project_id} for minimum level {min_level}")
            if min_level not in permission_levels:
                raise ValueError(f"Invalid minimum permission level {min_level} for project {project_id}")
            try:
                project = dxpy.DXProject(project_id)
                desc = await asyncio.to_thread(project.describe)
                user_perms = desc.get("level")
                if user_perms is None or user_perms == "NONE":
                    raise PermissionError(f"User '{user}' has no access to project '{project_id}'")
                if permission_levels[user_perms] < permission_levels[min_level]:
                    raise PermissionError(
                        f"User '{user}' does not have sufficient access to project '{project_id}'"
                        f" (has {user_perms}, needs {min_level})"
                    )
                logger.info(f"User '{user}' has '{user_perms}' access to project '{project_id}'")
                _DXPY_CONFIGURED[(project_id, min_level)] = True
            except dxpy.exceptions.DXError as e:
                logger.error(f"Could not access project '{project_id}': {e}")
                raise


async def submit_dx_job(
    applet_id: str,
    applet_project_id: str,
    project_id: str,
    input_dict: dict,
    instance_type: str,
    priority: str = "high",
    tags: list[str] | None = None,
    delay_seconds: int = 0,
) -> str:
    """Submit an applet job.

    Args:
        applet_id: Applet ID.
        applet_project_id: Project holding the applet.
        project_id: Project to run the job in.
        input_dict: Applet inputs, with files as dxlinks.
        instance_type: Instance type.
        priority: Job priority: `"high"`, `"normal"` or `"low"`.
        tags: Job tags.
        delay_seconds: Seconds to wait before submitting.

    Returns:
        The job ID.
    """
    if delay_seconds:
        await asyncio.sleep(delay_seconds)
    applet = dxpy.DXApplet(applet_id, project=applet_project_id)
    logger.debug(f"Submitting dx job {tags} ({applet_id}) with inputs {input_dict}")
    async with _DXAPI_SEMAPHORE:
        job = await asyncio.to_thread(
            applet.run,
            applet_input=input_dict,
            project=project_id,
            tags=tags,
            instance_type=instance_type,
            priority=priority,
        )
    return job.get_id()


async def wait_for_dx_job(
    job_id: str,
    poll_interval: float = 180,
    max_errors: int = 3,
    description: str = "",
    initial_delay: float | None = None,
) -> tuple[str, bool]:
    """Poll a job until it is done, failed or terminated.

    A job waiting to start is polled at doubling intervals, up to 16 times `poll_interval`.

    Args:
        job_id: Job ID.
        poll_interval: Seconds between polls.
        max_errors: Consecutive polling errors tolerated before giving up.
        description: Job name, for logs.
        initial_delay: Seconds before the first poll; None for one `poll_interval`.

    Returns:
        `(final_state, success)`; the state is `"polling_error"` after too many errors.
    """
    logger.debug(f"Polling dx job {description!r} ({job_id}), interval {poll_interval}s")
    await asyncio.sleep(poll_interval if initial_delay is None else initial_delay)
    job = dxpy.DXJob(job_id)
    consecutive_errors = 0
    consecutive_runnable = 0
    while True:
        try:
            async with _DXAPI_SEMAPHORE:
                describe = await asyncio.to_thread(job.describe)
            state = describe["state"]
            if state in ["done", "failed", "terminated"]:
                return state, state == "done"
            if state == "runnable":
                backoff = min(poll_interval * (2**consecutive_runnable), poll_interval * 16)
                consecutive_runnable += 1
                logger.debug(
                    f"Job {description!r} ({job_id}) is waiting to start, not polling for additional {backoff}s"
                )
                await asyncio.sleep(backoff)
            else:
                consecutive_runnable = 0
                logger.debug(f"Job {description!r} ({job_id}) in state '{state}', waiting…")
            consecutive_errors = 0
        except dxpy.exceptions.DXError as e:
            consecutive_errors += 1
            if consecutive_errors > max_errors:
                logger.error(
                    f"{consecutive_errors} consecutive errors polling job {description!r} ({job_id}): {e}"
                )
                return "polling_error", False
        await asyncio.sleep(poll_interval)


async def terminate_dx_job(job_id: str) -> None:
    """Terminate a job.

    Args:
        job_id: Job ID.
    """
    # Not _DXAPI_SEMAPHORE: an abort covers every running job at once and should not
    # queue behind polling; dxpy retries rate-limited calls
    async with _DXAPI_TERMINATE_SEMAPHORE:
        await asyncio.to_thread(dxpy.DXJob(job_id).terminate)


async def describe_jobs(job_ids: list[str], fields: dict[str, bool]) -> list[dict]:
    """Describe jobs with `/system/describeExecutions`, in batches.

    Args:
        job_ids: Job IDs.
        fields: Fields to describe, as `{name: True}`.

    Returns:
        The descriptions, omitting jobs the API cannot describe.
    """
    descriptions: list[dict] = []
    for i in range(0, len(job_ids), _DESCRIBE_BATCH_SIZE):
        async with _DXAPI_SEMAPHORE:
            response = await asyncio.to_thread(
                dxpy.api.system_describe_executions,
                {"executions": job_ids[i : i + _DESCRIBE_BATCH_SIZE], "fields": fields},
            )
        descriptions.extend(r["describe"] for r in response.get("results", []) if "describe" in r)
    return descriptions


def upload_file(
    file: Path,
    destination_folder: str,
    project_id: str,
    overwrite: bool = False,
) -> dxpy.bindings.dxfile.DXFile:
    """Upload a file to a project folder, creating the folder if needed.

    Args:
        file: Local file.
        destination_folder: Project folder.
        project_id: Project ID.
        overwrite: Remove files of the same name in the folder first.

    Returns:
        The uploaded file's handler.
    """
    project = dxpy.DXProject(project_id)
    project.new_folder(destination_folder, parents=True)

    if overwrite:
        existing = dxpy.find_data_objects(
            name=file.name,
            folder=destination_folder,
            project=project_id,
            classname="file",
            return_handler=True,
            recurse=False,
        )
        for f in existing:
            f.remove()

    return dxpy.upload_local_file(
        filename=str(file.absolute()),
        project=project_id,
        folder=destination_folder,
        wait_on_close=True,
        multithread=False,
    )


async def upload_file_async(
    file: Path,
    destination_folder: str,
    project_id: str,
    overwrite: bool = False,
) -> dxpy.bindings.dxfile.DXFile:
    """Run [`upload_file`][fedcat.platforms.dnanexus.dx.upload_file] in a thread."""
    async with _DXAPI_SEMAPHORE:
        return await asyncio.to_thread(upload_file, file, destination_folder, project_id, overwrite)


def download_files(
    project: str,
    folder: str,
    pattern: str,
    destination_folder: Path,
    max_workers: int | None = None,
) -> None:
    """Download the files of a project folder, and its subfolders, whose names match a regex.

    Args:
        project: Project ID.
        folder: Project folder.
        pattern: Regex the file names must match.
        destination_folder: Local directory, created if needed.
        max_workers: Most concurrent downloads, capped at 16; None for `NCPU * 4`.

    Raises:
        RuntimeError: If several matching files share a name, or any download fails.
    """
    logger.debug(f"Downloading from {project}:{folder} matching {pattern!r} → {destination_folder}")
    destination_folder.mkdir(parents=True, exist_ok=True)
    files = list(
        dxpy.find_data_objects(
            classname="file",
            project=project,
            folder=folder,
            name=pattern,
            name_mode="regexp",
            return_handler=True,
        )
    )
    # DNAnexus allows repeated names; they would land on the same local path, stale or not
    repeated = sorted(name for name, n in Counter(f.name for f in files).items() if n > 1)
    if repeated:
        raise RuntimeError(
            f"Several files share each of these names in {project}:{folder}, remove the stale copies: {repeated}"
        )

    if max_workers is None:
        max_workers = NCPU * 4
    max_workers = min(max_workers, 16)

    def _download_one(dxfile):
        """Download one file, returning `(name, error or None)`."""
        local_path = destination_folder / dxfile.name
        try:
            dxpy.download_dxfile(dxfile, str(local_path), project=project)
            return dxfile.name, None
        except Exception as e:  # noqa: BLE001 — returned and reported with the other failed downloads
            return dxfile.name, e

    failures: list[tuple[str, Exception]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as exe:
        futures = {exe.submit(_download_one, f): f for f in files}
        for fut in concurrent.futures.as_completed(futures):
            name, err = fut.result()
            if err:
                logger.error(f"Failed to download {name}: {err}")
                failures.append((name, err))
            else:
                logger.debug(f"Downloaded {name}")

    if failures:
        raise RuntimeError(f"Failed to download {len(failures)} files: {[n for n, _ in failures]}")


async def download_files_async(
    project: str,
    folder: str,
    pattern: str,
    destination_folder: Path,
    max_workers: int | None = None,
) -> None:
    """Run [`download_files`][fedcat.platforms.dnanexus.dx.download_files] in a thread, holding CPU permits."""
    permits = math.ceil(max_workers / 4) if max_workers is not None else None
    async with Permits(permits):
        await asyncio.to_thread(download_files, project, folder, pattern, destination_folder, max_workers)


def list_files(
    project: str,
    folder: str,
    pattern: str | None = None,
) -> list[str]:
    """List the file names in a project folder and its subfolders.

    Args:
        project: Project ID.
        folder: Project folder.
        pattern: Regex the names must match from their start; None for all files.

    Returns:
        Unique file names, sorted.
    """
    logger.debug(f"Listing files from {project}:{folder} with pattern {pattern!r}")
    files = dxpy.find_data_objects(
        classname="file",
        project=project,
        folder=folder,
        return_handler=True,
    )
    regex = re.compile(pattern) if pattern else None
    out: set[str] = set()
    for dxfile in files:
        name = dxfile.name
        if regex is not None and not regex.match(name):
            continue
        out.add(name)
    return sorted(out)


async def list_files_async(
    project: str,
    folder: str,
    pattern: str | None = None,
) -> list[str]:
    """Run [`list_files`][fedcat.platforms.dnanexus.dx.list_files] in a thread."""
    async with _DXAPI_SEMAPHORE:
        return await asyncio.to_thread(list_files, project, folder, pattern)


def file_id_to_path(file_id: str, project: str) -> str:
    """Resolve a file ID to its path.

    Args:
        file_id: `file-xxxx` ID.
        project: Project ID.

    Returns:
        The `/folder/name` path.

    Raises:
        dxpy.exceptions.DXError: If the file does not exist in the project.
    """
    desc: dict = dxpy.describe(file_id, fields={"folder": True, "name": True}, project=project)  # type: ignore
    folder = desc["folder"].rstrip("/")
    return f"{folder}/{desc['name']}"


def path_to_file_id(path: str, project: str) -> str:
    """Resolve a path to its file ID.

    Args:
        path: `/folder/name` path.
        project: Project ID.

    Returns:
        The `file-xxxx` ID.

    Raises:
        dxpy.exceptions.DXSearchError: If no file, or more than one, has the path.
    """
    parts = path.rsplit("/", 1)
    if len(parts) == 2:
        folder, name = parts
        folder = folder or "/"
    else:
        folder, name = "/", parts[0]
    try:
        result: dict = dxpy.find_one_data_object(
            classname="file",
            name=name,
            folder=folder,
            project=project,
            zero_ok=False,
            more_ok=False,
        )  # type: ignore
        return result["id"]
    except (dxpy.DXError, dxpy.exceptions.DXSearchError) as e:
        logger.error(f"Error resolving path {path!r} to file ID in project {project}: {e}")
        raise


def get_file_path(file_id: str, project: str) -> str:
    """Return the `/folder/name` path of a file given by ID or path.

    Args:
        file_id: `file-xxxx` ID or `/folder/name` path.
        project: Project ID.
    """
    if file_id.startswith("file-"):
        return file_id_to_path(file_id, project)
    else:
        return file_id


def get_file_id(path_or_id: str, project: str) -> str:
    """Return the `file-xxxx` ID of a file given by ID or path.

    Args:
        path_or_id: `file-xxxx` ID or `/folder/name` path.
        project: Project ID.
    """
    if path_or_id.startswith("file-"):
        return path_or_id
    else:
        return path_to_file_id(path_or_id, project)


def count_files_at_path(path: str, project: str) -> int:
    """Count the files with a path, up to 2.

    Args:
        path: `/folder/name` path.
        project: Project ID.

    Returns:
        0, 1, or 2 for more than one.
    """
    folder, _, name = path.rpartition("/")
    return sum(
        1
        for _ in dxpy.find_data_objects(
            classname="file", name=name, folder=folder or "/", project=project, recurse=False, limit=2
        )
    )


def file_exists(path_or_id: str, project: str) -> bool:
    """Return whether a file exists, and for a path, that it is unique.

    Args:
        path_or_id: `file-xxxx` ID or `/folder/name` path.
        project: Project ID.
    """
    try:
        if path_or_id.startswith("file-"):
            dxpy.describe(path_or_id, fields={"id": True}, project=project)
        else:
            path_to_file_id(path_or_id, project)
        return True
    except (dxpy.exceptions.DXError, dxpy.exceptions.DXSearchError, ValueError):
        return False
