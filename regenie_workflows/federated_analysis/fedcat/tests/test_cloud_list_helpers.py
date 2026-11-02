"""Unit tests for cloud file-list helper functions."""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
from fedcat.platforms.dnanexus import dx
from fedcat.platforms.gcs_cromwell import gcs


class _DummyDXFile:
    """Stand-in for a dxpy file handler."""

    def __init__(self, name: str):
        self.name = name


class _DummyBlob:
    """Stand-in for a GCS blob."""

    def __init__(self, name: str):
        self.name = name


class _DummyBucket:
    """Stand-in for a GCS bucket listing fixed blob names."""

    def __init__(self, names: list[str]):
        self._names = names

    def list_blobs(self, prefix: str):
        return [_DummyBlob(name) for name in self._names if name.startswith(prefix)]


class _DummyClient:
    """Stand-in for a GCS client whose buckets list fixed blob names."""

    def __init__(self, names: list[str]):
        self._names = names

    def bucket(self, _bucket_name: str):
        return _DummyBucket(self._names)


def test_list_ukb_files_filters_with_regex(monkeypatch):
    names = ["a.output_filelist", "b.regenie.gz", "c.output_filelist", "c.output_filelist"]

    def _fake_find_data_objects(**_kwargs):
        return [_DummyDXFile(name) for name in names]

    monkeypatch.setattr("fedcat.platforms.dnanexus.dx.dxpy.find_data_objects", _fake_find_data_objects)

    got = dx.list_files(
        project="project-123",
        folder="/tmp/out",
        pattern=r".*\.output_filelist$",
    )

    assert got == ["a.output_filelist", "c.output_filelist"]


def test_list_aou_files_filters_with_regex(monkeypatch):
    names = [
        "prefix/a.output_filelist",
        "prefix/b.regenie.gz",
        "prefix/c.output_filelist",
        "other/d.output_filelist",
    ]

    monkeypatch.setattr("fedcat.platforms.gcs_cromwell.gcs.storage.Client", lambda: _DummyClient(names))

    got = gcs.list_files(folder="bucket/prefix", pattern=r".*\.output_filelist$")

    assert got == ["a.output_filelist", "c.output_filelist"]


def test_list_aou_files_no_pattern_returns_all_basenames(monkeypatch):
    names = [
        "prefix/a.output_filelist",
        "prefix/b.regenie.gz",
        "prefix/c.output_filelist",
    ]

    monkeypatch.setattr("fedcat.platforms.gcs_cromwell.gcs.storage.Client", lambda: _DummyClient(names))

    got = gcs.list_files(folder="bucket/prefix")

    assert got == ["a.output_filelist", "b.regenie.gz", "c.output_filelist"]


def test_list_aou_files_ignores_sibling_folders_sharing_the_prefix(monkeypatch):
    names = [
        "run/agg_cond_step1/a.workflow.log",
        "run/agg_cond_step10/b.workflow.log",
        "run/agg_p0.0001/c.log",
        "run/agg_p0.00015/d.log",
    ]

    monkeypatch.setattr("fedcat.platforms.gcs_cromwell.gcs.storage.Client", lambda: _DummyClient(names))

    assert gcs.list_files(folder="bucket/run/agg_cond_step1") == ["a.workflow.log"]
    assert gcs.list_files(folder="bucket/run/agg_p0.0001/") == ["c.log"]


def test_dx_download_refuses_repeated_names(monkeypatch, tmp_path):
    """Two files with one name would race to the same local path, so the download stops instead."""
    import pytest

    names = ["a.tar.gz", "a.tar.gz", "b.tar.gz"]
    monkeypatch.setattr(
        "fedcat.platforms.dnanexus.dx.dxpy.find_data_objects", lambda **_kw: [_DummyDXFile(n) for n in names]
    )
    with pytest.raises(RuntimeError, match=r"\['a.tar.gz'\]"):
        dx.download_files(project="project-1", folder="/out", pattern=".*", destination_folder=tmp_path)


def test_dx_exists_reports_repeated_names(monkeypatch):
    import asyncio

    import pytest

    from fedcat.platforms.base import TransferError
    from fedcat.platforms.dnanexus.platform import DxPlatform

    found: list[str] = []
    monkeypatch.setattr(
        "fedcat.platforms.dnanexus.dx.dxpy.find_data_objects",
        lambda **kw: [_DummyDXFile(kw["name"])] * len(found),
    )
    platform = DxPlatform.__new__(DxPlatform)
    platform._cfg = type("Cfg", (), {"project_id": "project-1"})()
    for n, expected in ((0, False), (1, True)):
        found[:] = ["x"] * n
        assert asyncio.run(platform.exists("/out/a.tsv")) is expected
    found[:] = ["x", "x"]
    with pytest.raises(TransferError, match="several files"):
        asyncio.run(platform.exists("/out/a.tsv"))
