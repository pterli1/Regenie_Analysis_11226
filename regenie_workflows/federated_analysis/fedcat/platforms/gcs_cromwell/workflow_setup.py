"""`fedcat setup-workflows`: stage the repo's WDLs and their dependencies into a Workbench workspace.

Defines which WDLs and scripts make up a biobank's workflow set, and stages them:

1. The scripts the WDL tasks copy at run time are uploaded to
   `gs://<bucket>/workflows/<set>/v<version>/dependencies/`.
2. Each WDL's `dependencies_path` default is set to that folder, and the WDL is
   uploaded to `.../v<version>/wdls/` and registered with `wb workflow create` as
   `<stem>-v<version>` (dots as dashes, e.g. `-v1-7-0`), so versions sit side by side.

Also finds the WDLs and workflow options a `local_cromwell` submission sends. Only AoU
is implemented. See the
[`wb workflow` reference](https://support.workbench.verily.com/docs/references/cli_reference/wb/workflow/).
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

from loguru import logger

from fedcat.core.structures import (
    DEFAULT_DEPENDENCIES_VERSION,
    WORKFLOW_SET_NAME,
    BioBank,
    compute_dependencies_path,
    versioned_workflow_name,
)

from .gcs import upload_file

# ── Workflow set manifest ─────────────────────────────────────────────────────
# Paths relative to the biobank_workflows repo root, found by _resolve_repo_root
_AOU_WDL_DIR_REL = "regenie_workflows/aou_regenie_workflows/wdls"
_AOU_DEPENDENCIES = [
    # Run by the WDL tasks
    "regenie_workflows/aou_regenie_workflows/scripts/regenie_conditioning_agg_aou.sh",
    "regenie_workflows/aou_regenie_workflows/scripts/regenie_conditioning_singles_aou.sh",
    "regenie_workflows/aou_regenie_workflows/scripts/regenie_step2_aggregates_aou.sh",
    "regenie_workflows/aou_regenie_workflows/scripts/regenie_step2_single_variants_aou.sh",
    # Run by those scripts
    "regenie_workflows/shared_resources/run_regenie_step2_aggregates.sh",
    "regenie_workflows/shared_resources/run_regenie_step2_single_variants.sh",
    # Sourced helper functions
    "regenie_workflows/shared_resources/shared_functions.sh",
    "regenie_workflows/shared_resources/annotation_helpers.sh",
    # cond_engine "mancojo"
    "regenie_workflows/shared_resources/run_mancojo_conditional.sh",
]

# Staged when present, warned about when not: the manc-cojo binary is git-ignored, and
# only needed for cond_engine "mancojo"
_AOU_OPTIONAL_DEPENDENCIES = [
    "regenie_workflows/shared_resources/manc_cojo",
]

# Directories identifying the repo root
_REPO_MARKERS = (_AOU_WDL_DIR_REL, "regenie_workflows/shared_resources")


def _is_repo_root(candidate: Path) -> bool:
    """Return whether a directory is the biobank_workflows repo root."""
    return all((candidate / marker).is_dir() for marker in _REPO_MARKERS)


def _auto_detect_repo_root() -> Path | None:
    """Return the biobank_workflows repo root above the fedcat package or working directory, or None."""
    seen: set[Path] = set()
    for start in (Path(__file__).resolve(), Path.cwd().resolve()):
        for candidate in (start, *start.parents):
            if candidate in seen:
                continue
            seen.add(candidate)
            if _is_repo_root(candidate):
                return candidate
    return None


def _resolve_repo_root(explicit: str | Path | None = None, *, interactive: bool = True) -> Path:
    """Return the biobank_workflows repo root: the given path, else found, else asked for.

    Raises:
        WorkflowSetupError: If the given path is not the repo root, or none is found.
    """
    if explicit:
        root = Path(explicit).expanduser().resolve()
        if not _is_repo_root(root):
            raise WorkflowSetupError(
                f"--repo-root {root} does not look like the biobank_workflows repo "
                f"(missing {' / '.join(_REPO_MARKERS)})."
            )
        return root

    detected = _auto_detect_repo_root()
    if detected is not None:
        return detected

    if interactive and sys.stdin.isatty():
        while True:
            entered = input(
                "Could not locate the biobank_workflows repo — enter its path (blank to abort): "
            ).strip()
            if not entered:
                break
            candidate = Path(entered).expanduser().resolve()
            if _is_repo_root(candidate):
                return candidate
            print(f"  {candidate} is missing {' / '.join(_REPO_MARKERS)}; try again.")

    raise WorkflowSetupError(
        "Could not locate the biobank_workflows repo (searched from the fedcat package "
        "and the current directory). Pass --repo-root explicitly."
    )


# ── Local Cromwell submission sources ─────────────────────────────────────────
# A local Cromwell server is sent the WDL and workflow options with every submission.
# Copies in the working directory take priority, to change them without editing the repo
LOCAL_OPTIONS_NAME = "cromwell_options.json"
DEFAULT_WORKFLOW_OPTIONS = {
    "workflow_failure_mode": "ContinueWhilePossible",
    "write_to_cache": False,
    "read_from_cache": False,
}
_ANNOUNCED_LOCAL: set[str] = set()


def _announce_local(name: str) -> None:
    """Log, once per file, that a working-directory copy is used."""
    if name not in _ANNOUNCED_LOCAL:
        _ANNOUNCED_LOCAL.add(name)
        logger.info(f"Using {name} from the working directory {Path.cwd()}")


def resolve_local_wdl(stem: str) -> Path | None:
    """Return the WDL to submit to a local server: the working directory's copy, else the repo's.

    Args:
        stem: WDL workflow name.

    Returns:
        The WDL path, or None if neither exists.
    """
    local = Path(f"{stem}.wdl")
    if local.is_file():
        _announce_local(local.name)
        return local
    root = _auto_detect_repo_root()
    if root is not None and (bundled := root / _AOU_WDL_DIR_REL / local.name).is_file():
        return bundled
    return None


def local_workflow_options() -> str:
    """Return the workflow options JSON: the working directory's `cromwell_options.json`, else `DEFAULT_WORKFLOW_OPTIONS`."""
    local = Path(LOCAL_OPTIONS_NAME)
    if local.is_file():
        _announce_local(local.name)
        return local.read_text()
    return json.dumps(DEFAULT_WORKFLOW_OPTIONS)


# The WDLs' dependencies_path = "gs://..." defaults; group 1 is up to the opening quote
_DEPENDENCIES_DEFAULT_RE = re.compile(
    r'(dependencies_path\s*=\s*")gs://[^"]*"',
)


class WorkflowSetupError(RuntimeError):
    """Raised when staging or registering a workflow set fails."""


def _require_aou(biobank: BioBank) -> None:
    """Raise NotImplementedError for biobanks other than AoU."""
    if biobank is not BioBank.AOU:
        raise NotImplementedError(f"Workflow staging is only implemented for AoU, not {biobank.value!r}")


def _dependency_paths(biobank: BioBank, repo_root: Path) -> list[Path]:
    """Return the dependency files to stage, raising WorkflowSetupError if a required one is missing."""
    _require_aou(biobank)
    paths = [repo_root / rel for rel in _AOU_DEPENDENCIES]
    missing = [str(p) for p in paths if not p.is_file()]
    if missing:
        raise WorkflowSetupError("Dependency file(s) not found:\n  " + "\n  ".join(missing))
    for rel in _AOU_OPTIONAL_DEPENDENCIES:
        optional = repo_root / rel
        if optional.is_file():
            paths.append(optional)
        else:
            logger.warning(
                f"Optional dependency {rel} not found; staging without it. "
                "cond_engine='mancojo' runs will fail until it is added."
            )
    return paths


def _wdl_paths(biobank: BioBank, repo_root: Path) -> list[Path]:
    """Return the WDLs to stage, raising WorkflowSetupError if there are none."""
    _require_aou(biobank)
    wdl_dir = repo_root / _AOU_WDL_DIR_REL
    wdls = sorted(wdl_dir.glob("*.wdl"))
    if not wdls:
        raise WorkflowSetupError(f"No .wdl files found in {wdl_dir}")
    return wdls


def dependency_filenames(biobank: BioBank, *, cond_engine: str = "regenie") -> list[str]:
    """Return the names of the files a run needs in the staged dependencies folder.

    Taken from the manifest, so needs no repo checkout.

    Args:
        biobank: Biobank; only AoU is supported.
        cond_engine: The run's `cond_engine`; `"mancojo"` adds the manc-cojo binary.

    Returns:
        File names.

    Raises:
        NotImplementedError: If `biobank` is not AoU.
    """
    _require_aou(biobank)
    rels = list(_AOU_DEPENDENCIES)
    if cond_engine == "mancojo":
        rels += _AOU_OPTIONAL_DEPENDENCIES
    return [Path(rel).name for rel in rels]


def _rewrite_dependencies_default(wdl_text: str, new_path: str) -> str:
    """Set every `dependencies_path` default in a WDL's text to `new_path`."""
    return _DEPENDENCIES_DEFAULT_RE.sub(rf'\g<1>{new_path}"', wdl_text)


def _run_wb(*args: str) -> subprocess.CompletedProcess[str]:
    """Run a `wb` command, raising WorkflowSetupError with `wb`'s message if it fails."""
    try:
        return subprocess.run(["wb", *args], check=True, text=True, capture_output=True)
    except subprocess.CalledProcessError as exc:
        detail = (exc.stderr or exc.stdout or "").strip()
        raise WorkflowSetupError(f"`wb {' '.join(args)}` failed (exit {exc.returncode}): {detail}") from exc


def registered_workflows() -> list[dict]:
    """Return the workflows registered in the workspace, as `wb workflow list` describes them.

    Raises:
        WorkflowSetupError: If the workflows cannot be listed.
    """
    try:
        result = _run_wb("workflow", "list", "--format=JSON")
        workflows = json.loads(result.stdout)
    except (WorkflowSetupError, json.JSONDecodeError, FileNotFoundError) as exc:
        raise WorkflowSetupError(f"Failed to list existing workflows via 'wb workflow list': {exc}") from exc
    return [wf for wf in workflows if isinstance(wf, dict) and wf.get("id")]


def wdl_bucket_path(wdl_name: str, version: str) -> str:
    """Return where `setup-workflows` stages a WDL at a version, relative to the workspace bucket."""
    return f"workflows/{WORKFLOW_SET_NAME}/v{version}/wdls/{wdl_name}.wdl"


def workspace_bucket_uuid(bucket: str, bucket_resource_id: str = "") -> str | None:
    """Return the workspace bucket's UUID, by which registered workflows refer to it.

    Args:
        bucket: Bucket name, with or without `gs://`.
        bucket_resource_id: Its Workbench resource id, also matched.

    Returns:
        The UUID, or None if it cannot be found.
    """
    try:
        resources = json.loads(_run_wb("resource", "list", "--type=GCS_BUCKET", "--format=JSON").stdout)
    except (WorkflowSetupError, json.JSONDecodeError, FileNotFoundError):
        return None
    name = bucket.removeprefix("gs://").strip("/")
    for res in resources:
        if isinstance(res, dict) and (res.get("bucketName") == name or res.get("id") == bucket_resource_id):
            return res.get("uuid") or None
    return None


def resolve_registered_workflow(
    wdl_name: str, version: str, registered: list[dict], bucket_uuid: str | None = None
) -> str:
    """Return the id of the registered workflow running a WDL as staged at a version.

    Matched on the staged WDL's path, so the workflow may have any name. Users staging
    the same version to their own buckets share the path, so only workflows from this
    bucket count when its UUID is known.

    Args:
        wdl_name: WDL workflow name.
        version: Dependencies version.
        registered: From [`registered_workflows`][fedcat.platforms.gcs_cromwell.workflow_setup.registered_workflows].
        bucket_uuid: Workspace bucket UUID, or None if unknown.

    Returns:
        The workflow id.

    Raises:
        WorkflowSetupError: If no workflow, or more than one, runs the WDL.
    """
    path = wdl_bucket_path(wdl_name, version)
    matches = [wf for wf in registered if (wf.get("bucketSource") or {}).get("path") == path]
    if bucket_uuid:
        matches = [wf for wf in matches if wf["bucketSource"].get("bucketResourceId") == bucket_uuid]
    if len(matches) == 1:
        return matches[0]["id"]
    if not matches:
        raise WorkflowSetupError(
            f"no workflow registered for {path} in this workspace's bucket; "
            f"run 'fedcat setup-workflows --version {version}'"
        )
    found = ", ".join(f"{wf['id']} (by {wf.get('createdBy') or '?'})" for wf in matches)
    raise WorkflowSetupError(
        f"several workflows run {path} and the workspace bucket did not tell them apart: {found}"
    )


def _confirm_replace(workflow_id: str, existing: dict) -> bool:
    """Ask whether to replace an existing workflow."""
    source = existing.get("bucketSource") or {}
    answer = input(
        f"    '{workflow_id}' already exists (registered by {existing.get('createdBy') or '?'} "
        f"for {source.get('path') or '?'}). Replace it? [y/N] "
    )
    return answer.strip().lower() in ("y", "yes")


def setup_workflows(
    *,
    biobank: BioBank = BioBank.AOU,
    bucket: str,
    bucket_id: str,
    version: str = DEFAULT_DEPENDENCIES_VERSION,
    repo_root: str | Path | None = None,
    force: bool = False,
    skip_existing: bool = False,
    dry_run: bool = False,
) -> None:
    """Stage a biobank's workflow set into the workspace and register its workflows.

    Args:
        biobank: Biobank; only AoU is supported.
        bucket: Workspace bucket name, with or without `gs://`, to upload to.
        bucket_id: Its Workbench resource id, for `wb workflow create`.
        version: Dependencies version to stage as.
        repo_root: biobank_workflows checkout; found, or asked for, when None.
        force: Replace existing workflows without asking.
        skip_existing: Keep existing workflows without asking.
        dry_run: Log what would be done without doing it.

    Raises:
        WorkflowSetupError: If both `force` and `skip_existing` are given, the repo or
            files are not found, or a `wb` command fails.
        NotImplementedError: If `biobank` is not AoU.
    """
    _require_aou(biobank)
    if force and skip_existing:
        raise WorkflowSetupError("force and skip_existing are mutually exclusive")

    root = _resolve_repo_root(repo_root)
    bucket_name = bucket.removeprefix("gs://").strip("/")
    dep_path = compute_dependencies_path(bucket, version)
    dep_files = _dependency_paths(biobank, root)
    wdls = _wdl_paths(biobank, root)

    logger.info(f"Staging '{WORKFLOW_SET_NAME}' v{version} for {biobank.value} from {root} -> gs://{bucket_name}")

    # 1. Upload the dependencies
    dep_prefix = f"{bucket_name}/workflows/{WORKFLOW_SET_NAME}/v{version}/dependencies"
    logger.info(f"Uploading {len(dep_files)} dependency file(s) -> gs://{dep_prefix}/")
    for dep in dep_files:
        dest = f"{dep_prefix}/{dep.name}"
        if dry_run:
            logger.info(f"    [dry-run] upload {dep.name} -> gs://{dest}")
        else:
            upload_file(dep, dest)
            logger.debug(f"    uploaded {dep.name}")

    # 2. Upload each WDL, with its dependencies_path set, and register it
    existing = {} if dry_run else {wf["id"]: wf for wf in registered_workflows()}
    logger.info(f"Registering {len(wdls)} workflow(s) from {root / _AOU_WDL_DIR_REL}")
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        for wdl in wdls:
            # Versioned, so staging a version never replaces one someone else is running
            workflow_id = versioned_workflow_name(wdl.stem, version)
            logger.info(f"- {workflow_id}")

            if workflow_id in existing:
                if skip_existing:
                    logger.info("    already present -> skipping (skip_existing)")
                    continue
                if not (force or _confirm_replace(workflow_id, existing[workflow_id])):
                    logger.info("    keeping existing workflow -> skipping")
                    continue
                if dry_run:
                    logger.info(f"    [dry-run] wb workflow delete --workflow={workflow_id}")
                else:
                    _run_wb("workflow", "delete", f"--workflow={workflow_id}", "--quiet")

            baked = _rewrite_dependencies_default(wdl.read_text(), dep_path)
            if "dependencies_path" in baked and dep_path not in baked:
                logger.warning(f"    {wdl.name}: dependencies_path default not rewritten (unexpected format)")
            bucket_path = wdl_bucket_path(wdl.stem, version)

            if dry_run:
                logger.info(f"    [dry-run] upload (baked) {wdl.name} -> gs://{bucket_name}/{bucket_path}")
                logger.info(
                    f"    [dry-run] wb workflow create --workflow={workflow_id} "
                    f"--bucket-id={bucket_id} --path={bucket_path}"
                )
                continue

            staged = tmp_dir / wdl.name
            staged.write_text(baked)
            upload_file(staged, f"{bucket_name}/{bucket_path}")
            _run_wb(
                "workflow",
                "create",
                f"--workflow={workflow_id}",
                f"--display-name={workflow_id}",
                f"--bucket-id={bucket_id}",
                f"--path={bucket_path}",
            )
            logger.debug(f"    registered {workflow_id}")

    logger.info("Workflow staging complete.")
