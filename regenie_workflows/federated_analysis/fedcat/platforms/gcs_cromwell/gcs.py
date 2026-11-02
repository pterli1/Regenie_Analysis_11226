"""Google Cloud Storage helpers, with paths given as `bucket/path` without `gs://`."""

from __future__ import annotations

import asyncio
import re
from pathlib import Path

from google.cloud import storage

from fedcat.core.logging import Permits


def upload_file(file: Path, destination: str) -> None:
    """Upload a file to GCS.

    Args:
        file: Local file.
        destination: Object path, `bucket/path/to/object`.
    """
    client = storage.Client()
    parts = destination.split("/")
    bucket = client.bucket(parts[0])
    blob = bucket.blob("/".join(parts[1:]))
    blob.upload_from_filename(str(file.absolute()))


async def upload_file_async(file: Path, destination: str) -> None:
    """Run [`upload_file`][fedcat.platforms.gcs_cromwell.gcs.upload_file] in a thread."""
    await asyncio.to_thread(upload_file, file, destination)


def _bucket_and_prefix(folder: str) -> tuple[str, str]:
    """Split `bucket/path` into the bucket and an object prefix ending in `/`.

    The trailing `/` stops a folder matching others it is a prefix of (`step1` and `step10`).
    """
    bucket_name, _, path = folder.strip("/").partition("/")
    return bucket_name, f"{path}/" if path else ""


def download_files(folder: str, pattern: str, destination_folder: Path) -> None:
    """Download the files under a GCS folder, at any depth, whose names match a regex.

    Args:
        folder: Folder, `bucket/path`.
        pattern: Regex the file names must match from their start.
        destination_folder: Local directory, created if needed.
    """
    destination_folder.mkdir(parents=True, exist_ok=True)
    client = storage.Client()
    bucket_name, prefix = _bucket_and_prefix(folder)
    bucket = client.bucket(bucket_name)
    blobs = [blob for blob in bucket.list_blobs(prefix=prefix) if re.match(pattern, blob.name.split("/")[-1])]
    for blob in blobs:
        blob.download_to_filename(destination_folder / blob.name.split("/")[-1])


async def download_files_async(folder: str, pattern: str, destination_folder: Path) -> None:
    """Run [`download_files`][fedcat.platforms.gcs_cromwell.gcs.download_files] in a thread."""
    async with Permits(1):
        await asyncio.to_thread(download_files, folder, pattern, destination_folder)


def list_files(folder: str, pattern: str | None = None) -> list[str]:
    """List the file names under a GCS folder, at any depth.

    Args:
        folder: Folder, `bucket/path`.
        pattern: Regex the names must match from their start; None for all files.

    Returns:
        Unique file names, sorted.
    """
    client = storage.Client()
    bucket_name, prefix = _bucket_and_prefix(folder)
    bucket = client.bucket(bucket_name)
    regex = re.compile(pattern) if pattern else None
    out: set[str] = set()
    for blob in bucket.list_blobs(prefix=prefix):
        name = blob.name.split("/")[-1]
        if regex is not None and not regex.match(name):
            continue
        out.add(name)
    return sorted(out)


async def list_files_async(folder: str, pattern: str | None = None) -> list[str]:
    """Run [`list_files`][fedcat.platforms.gcs_cromwell.gcs.list_files] in a thread."""
    async with Permits(1):
        return await asyncio.to_thread(list_files, folder, pattern)


def blob_exists(full_gcs_path: str) -> bool:
    """Return whether a GCS object exists.

    Args:
        full_gcs_path: Object path, `gs://bucket/path`.

    Raises:
        ValueError: If the path does not start with `gs://` or has no object path.
    """
    if not full_gcs_path.startswith("gs://"):
        raise ValueError(f"Expected gs:// path, got: {full_gcs_path!r}")
    stripped = full_gcs_path[len("gs://") :]
    parts = stripped.split("/", 1)
    if len(parts) != 2 or not parts[1]:
        raise ValueError(f"GCS path must include bucket and blob path: {full_gcs_path!r}")
    bucket_name, blob_path = parts
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    return bucket.blob(blob_path).exists()
