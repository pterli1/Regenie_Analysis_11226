"""Platforms: the storage and job services a cohort runs on.

- [`base`][fedcat.platforms.base]: the [`Platform`][fedcat.platforms.base.Platform]
  interface every platform implements.
- [`gcs_cromwell`][fedcat.platforms.gcs_cromwell]: GCS storage with Cromwell jobs, for AoU.
- [`dnanexus`][fedcat.platforms.dnanexus]: DNAnexus, for UKB.
- [`registry`][fedcat.platforms.registry]: builds each cohort's platform from
  `BIOBANK_PROFILES`.

A new platform is a new subpackage and a `BIOBANK_PROFILES` row; `fedcat.core` does not
depend on which biobank a cohort is.
"""

from fedcat.platforms.base import (
    AuthError,
    JobSpec,
    PathSyntaxError,
    Platform,
    PlatformError,
    SubmitError,
    TransferError,
)
from fedcat.platforms.registry import build_platforms

__all__ = [
    "AuthError",
    "JobSpec",
    "PathSyntaxError",
    "Platform",
    "PlatformError",
    "SubmitError",
    "TransferError",
    "build_platforms",
]
