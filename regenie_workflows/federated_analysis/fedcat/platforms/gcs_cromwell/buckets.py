"""Choosing the workspace bucket from the Verily Workbench workspace, with the `wb` CLI.

The CONTROLLED GCS buckets of the current workspace are listed with
`wb resource list --type=GCS_BUCKET` for the user to choose from. A bucket has two
identifiers:

- `resource_id`: the Workbench resource id (e.g. `anthropometric-traits`), passed to
  `wb` commands as `--bucket-id` or `--output-bucket-id`.
- `bucket_name`: the GCS bucket (e.g. `anthropometric-traits-wb-agile-asparagus-2878`),
  used in every `gs://` path.
"""

from __future__ import annotations

import json
import subprocess
import sys
from dataclasses import dataclass

from loguru import logger

_CONTROLLED = "CONTROLLED"


@dataclass(frozen=True)
class WorkspaceBucket:
    """A workspace bucket: its Workbench resource id and GCS bucket name."""

    resource_id: str
    """Workbench resource id."""
    bucket_name: str
    """GCS bucket name."""

    @property
    def gs_url(self) -> str:
        """`gs://<bucket_name>`."""
        return f"gs://{self.bucket_name}"


def list_controlled_buckets() -> list[WorkspaceBucket] | None:
    """List the CONTROLLED GCS buckets of the current workspace with `wb`.

    Returns:
        The buckets, or None if `wb` is unavailable or the listing fails.
    """
    try:
        result = subprocess.run(
            ["wb", "resource", "list", "--type=GCS_BUCKET", "--format=JSON"],
            check=True,
            text=True,
            capture_output=True,
        )
        resources = json.loads(result.stdout)
    except (subprocess.CalledProcessError, json.JSONDecodeError, FileNotFoundError, OSError) as exc:
        logger.debug(f"wb bucket listing unavailable: {exc}")
        return None

    buckets: list[WorkspaceBucket] = []
    for res in resources:
        if not isinstance(res, dict) or res.get("stewardshipType") != _CONTROLLED:
            continue
        bucket_name = res.get("bucketName")
        resource_id = res.get("id")
        if bucket_name and resource_id:
            buckets.append(WorkspaceBucket(resource_id=resource_id, bucket_name=bucket_name))
    return buckets


def resolve_configured_workspace_bucket(bucket: str, resource_id: str = "") -> WorkspaceBucket:
    """Find a configured bucket among the workspace's CONTROLLED buckets.

    Args:
        bucket: GCS bucket name, with or without `gs://`.
        resource_id: Its Workbench resource id, if configured; looked up when empty.

    Returns:
        The workspace bucket.

    Raises:
        ValueError: If `bucket` is empty or not a CONTROLLED bucket of the workspace,
            or `resource_id` is not its resource id.
        RuntimeError: If `wb` cannot list the buckets.
    """
    bucket_name = bucket.removeprefix("gs://").strip("/")
    if not bucket_name:
        raise ValueError("aou_bucket cannot be empty.")

    buckets = list_controlled_buckets()
    if buckets is None:
        raise RuntimeError(
            "Could not resolve the Workbench resource id for configured aou_bucket "
            f"{bucket!r}: `wb resource list` is unavailable."
        )

    matched = next((candidate for candidate in buckets if candidate.bucket_name == bucket_name), None)
    if matched is None:
        raise ValueError(f"Configured aou_bucket {bucket!r} is not a CONTROLLED GCS bucket in this workspace.")
    if resource_id and resource_id != matched.resource_id:
        raise ValueError(
            f"Configured aou_bucket {bucket!r} has resource id {matched.resource_id!r}, not {resource_id!r}."
        )
    return matched


def _prompt_manual_bucket() -> WorkspaceBucket:
    """Prompt for a bucket name and resource id, which defaults to the name."""
    bucket_name = input("Enter GCS bucket name (without gs://): ").strip().removeprefix("gs://").strip("/")
    if not bucket_name:
        raise ValueError("No bucket name entered.")
    resource_id = input(f"Enter the Workbench bucket resource id [default: {bucket_name}]: ").strip()
    return WorkspaceBucket(resource_id=resource_id or bucket_name, bucket_name=bucket_name)


def resolve_workspace_bucket(*, interactive: bool = True) -> WorkspaceBucket:
    """Ask the user to choose a workspace bucket from those `wb` lists, or to enter one.

    Without a terminal, the only bucket is chosen if there is just one.

    Args:
        interactive: Allow prompting; also needs stdin to be a terminal.

    Returns:
        The chosen bucket.

    Raises:
        RuntimeError: If there is no terminal and not exactly one bucket.
        ValueError: If an entered bucket name is empty.
    """
    interactive = interactive and sys.stdin.isatty()
    buckets = list_controlled_buckets()

    if not buckets:
        if buckets is None:
            logger.warning("Could not list workspace buckets via 'wb' (not installed or not authenticated).")
        else:
            logger.warning("No CONTROLLED GCS buckets found in the workspace.")
        if not interactive:
            raise RuntimeError(
                "No workspace bucket could be resolved and no interactive terminal is available. "
                "Pass the bucket explicitly (config aou_bucket / --bucket)."
            )
        return _prompt_manual_bucket()

    if len(buckets) == 1 and not interactive:
        logger.info(f"Auto-selecting the only CONTROLLED bucket: {buckets[0].resource_id}")
        return buckets[0]

    if not interactive:
        raise RuntimeError(
            f"Multiple workspace buckets found ({', '.join(b.resource_id for b in buckets)}) and no "
            "interactive terminal is available. Pass the bucket explicitly (config aou_bucket / --bucket)."
        )

    print("\nAvailable CONTROLLED workspace buckets:")
    for i, b in enumerate(buckets, start=1):
        print(f"  {i}. {b.resource_id}  (gs://{b.bucket_name})")
    print("  c. Enter a custom bucket")

    while True:
        choice = input(f"Select a bucket [1-{len(buckets)}, or c]: ").strip().lower()
        if choice == "c":
            return _prompt_manual_bucket()
        if choice.isdigit() and 1 <= int(choice) <= len(buckets):
            return buckets[int(choice) - 1]
        print("Invalid selection, try again.")
