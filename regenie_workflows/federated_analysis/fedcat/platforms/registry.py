"""Build each cohort's [`Platform`][fedcat.platforms.base.Platform], chosen by its `platform_kind`."""

from __future__ import annotations

from functools import cache

from fedcat.core.structures import BIOBANK_PROFILES, BioBank, DiscoveryConfig
from fedcat.platforms.base import Platform
from fedcat.platforms.dnanexus import DxCfg, DxPlatform
from fedcat.platforms.dnanexus.dx import UKB_APPLET_PROJECT
from fedcat.platforms.gcs_cromwell import GcsCromwellCfg, GcsCromwellPlatform


@cache
def _gcs_platform(cfg: GcsCromwellCfg) -> GcsCromwellPlatform:
    """Return the GCS/Cromwell platform for a config, one instance per distinct config."""
    return GcsCromwellPlatform(cfg)


@cache
def _dx_platform(cfg: DxCfg) -> DxPlatform:
    """Return the DNAnexus platform for a config, one instance per distinct config."""
    return DxPlatform(cfg)


def build_platforms(config: DiscoveryConfig) -> dict[BioBank, Platform]:
    """Build the platform of each cohort in a config.

    Cohorts with the same platform settings share one platform instance.

    Args:
        config: Run config.

    Returns:
        Platform per biobank.

    Raises:
        ValueError: If a biobank's `platform_kind` is unknown.
    """
    out: dict[BioBank, Platform] = {}
    for bb, cohort in config.cohorts.items():
        match BIOBANK_PROFILES[bb].platform_kind:
            case "gcs_cromwell":
                out[bb] = _gcs_platform(
                    GcsCromwellCfg(
                        bucket=config.aou_bucket,
                        bucket_resource_id=config.aou_bucket_resource_id,
                        backend=cohort.workflow_backend,
                        dependencies_version=cohort.dependencies_version,
                    )
                )
            case "dnanexus":
                # Cohort projects fall back to the top-level project and the shared applet project
                out[bb] = _dx_platform(
                    DxCfg(
                        project_id=cohort.dx_project_id or config.dx_project_id,
                        applet_project_id=cohort.applet_project_id or UKB_APPLET_PROJECT,
                    )
                )
            case other:
                raise ValueError(
                    f"Unknown platform_kind {other!r} for biobank {bb.value}; check BIOBANK_PROFILES"
                )
    return out
