"""Tests for resolving configured AoU workspace buckets."""

from __future__ import annotations

import pytest

from fedcat.platforms.gcs_cromwell.buckets import (
    WorkspaceBucket,
    resolve_configured_workspace_bucket,
)


def test_configured_bucket_derives_its_resource_id(monkeypatch):
    monkeypatch.setattr(
        "fedcat.platforms.gcs_cromwell.buckets.list_controlled_buckets",
        lambda: [WorkspaceBucket(resource_id="bucket-resource", bucket_name="workspace-bucket")],
    )

    resolved = resolve_configured_workspace_bucket("gs://workspace-bucket")

    assert resolved == WorkspaceBucket(resource_id="bucket-resource", bucket_name="workspace-bucket")


def test_configured_bucket_rejects_mismatched_resource_id(monkeypatch):
    monkeypatch.setattr(
        "fedcat.platforms.gcs_cromwell.buckets.list_controlled_buckets",
        lambda: [WorkspaceBucket(resource_id="bucket-resource", bucket_name="workspace-bucket")],
    )

    with pytest.raises(ValueError, match="not 'other-resource'"):
        resolve_configured_workspace_bucket("workspace-bucket", resource_id="other-resource")


def _wb_listing(monkeypatch, resources, returncode: int = 0):
    """Make `wb resource list` return `resources`, or fail with `returncode`."""
    import json
    import subprocess

    def _run(cmd, **_kwargs):
        if returncode:
            raise subprocess.CalledProcessError(returncode, cmd)
        return subprocess.CompletedProcess(cmd, 0, json.dumps(resources), "")

    monkeypatch.setattr("fedcat.platforms.gcs_cromwell.buckets.subprocess.run", _run)


# fields as `wb resource list --type=GCS_BUCKET --format=JSON` gives them
_MINE = {"uuid": "u1", "id": "Mine", "stewardshipType": "CONTROLLED", "bucketName": "mine-wb-x"}
_THEIRS = {"uuid": "u2", "id": "Theirs", "stewardshipType": "CONTROLLED", "bucketName": "theirs-wb-x"}
_REFERENCED = {"uuid": "u3", "id": "Ref", "stewardshipType": "REFERENCED", "bucketName": "ref"}


def test_listing_keeps_controlled_buckets(monkeypatch):
    from fedcat.platforms.gcs_cromwell.buckets import list_controlled_buckets

    _wb_listing(monkeypatch, [_MINE, _REFERENCED, "junk"])
    assert list_controlled_buckets() == [WorkspaceBucket(resource_id="Mine", bucket_name="mine-wb-x")]
    _wb_listing(monkeypatch, [], returncode=1)
    assert list_controlled_buckets() is None


def test_non_interactive_resolution_needs_a_single_bucket(monkeypatch):
    from fedcat.platforms.gcs_cromwell.buckets import resolve_workspace_bucket

    _wb_listing(monkeypatch, [_MINE])
    assert resolve_workspace_bucket(interactive=False).bucket_name == "mine-wb-x"
    _wb_listing(monkeypatch, [_MINE, _THEIRS])
    with pytest.raises(RuntimeError, match="Multiple workspace buckets"):
        resolve_workspace_bucket(interactive=False)
    _wb_listing(monkeypatch, [], returncode=1)
    with pytest.raises(RuntimeError, match="no interactive terminal"):
        resolve_workspace_bucket(interactive=False)
    # a configured bucket must be one of the workspace's
    _wb_listing(monkeypatch, [_MINE])
    with pytest.raises(ValueError, match="not a CONTROLLED GCS bucket"):
        resolve_configured_workspace_bucket("gs://elsewhere")
