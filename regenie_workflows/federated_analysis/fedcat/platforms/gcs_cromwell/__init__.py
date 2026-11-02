"""GCS and Cromwell platform, for All of Us and other Verily Workbench-hosted biobanks.

- [`gcs`][fedcat.platforms.gcs_cromwell.gcs]: GCS storage.
- [`cromwell`][fedcat.platforms.gcs_cromwell.cromwell],
  [`cromwell_client`][fedcat.platforms.gcs_cromwell.cromwell_client]: workflow
  submission and polling, through the Workbench CLI or a Cromwell server.
- [`cromwell_server`][fedcat.platforms.gcs_cromwell.cromwell_server]: `fedcat cromwell`.
- [`buckets`][fedcat.platforms.gcs_cromwell.buckets]: workspace bucket selection.
- [`workflow_setup`][fedcat.platforms.gcs_cromwell.workflow_setup]: `fedcat setup-workflows`.
- [`pricing`][fedcat.platforms.gcs_cromwell.pricing]: task cost estimates.
- [`GcsCromwellPlatform`][fedcat.platforms.gcs_cromwell.platform.GcsCromwellPlatform]:
  the platform interface.
"""

from fedcat.platforms.gcs_cromwell.platform import GcsCromwellCfg, GcsCromwellPlatform

__all__ = ["GcsCromwellCfg", "GcsCromwellPlatform"]
