"""Contract tests for the Platform abstraction.

Chiefly that `exists` raises `PathSyntaxError` for a malformed reference on every
platform, rather than reporting the file missing.
"""

from __future__ import annotations

import asyncio

import pytest

from fedcat.platforms.base import PathSyntaxError
from fedcat.platforms.dnanexus import DxCfg, DxPlatform
from fedcat.platforms.gcs_cromwell import GcsCromwellCfg, GcsCromwellPlatform


def _platforms():
    """Return one instance of each real platform."""
    return [
        GcsCromwellPlatform(GcsCromwellCfg(bucket="my-bucket")),
        DxPlatform(DxCfg(project_id="project-000000000000000000000000")),
    ]


@pytest.mark.parametrize("platform", _platforms(), ids=lambda p: p.kind)
def test_exists_raises_pathsyntaxerror_on_malformed_ref(platform):
    # A single garbage token is neither a valid GCS bucket/blob nor a DNAnexus
    # file-id / absolute path — both platforms must reject it as malformed rather
    # than reporting it as absent. The check happens before any network call.
    with pytest.raises(PathSyntaxError):
        asyncio.run(platform.exists("garbage"))


def test_dx_exists_rejects_bad_file_id_shape():
    platform = DxPlatform(DxCfg(project_id="project-000000000000000000000000"))
    # Looks file-ish but is not a 24-char id and not an absolute path.
    with pytest.raises(PathSyntaxError):
        asyncio.run(platform.exists("file-short"))


def test_platform_kinds_match_profiles():
    from fedcat.core.structures import BIOBANK_PROFILES

    kinds = {p.kind for p in _platforms()}
    profile_kinds = {prof.platform_kind for prof in BIOBANK_PROFILES.values()}
    assert profile_kinds <= kinds
