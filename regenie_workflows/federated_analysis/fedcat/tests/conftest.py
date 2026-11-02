"""Shared pytest fixtures and path setup for the fedcat test suite.

Puts the package root on `sys.path`, so tests import this checkout's `fedcat`.
"""

from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

import pytest

_REGENIE_WORKFLOWS = Path(__file__).resolve().parents[2]
if str(_REGENIE_WORKFLOWS) not in sys.path:
    sys.path.insert(0, str(_REGENIE_WORKFLOWS))

from fedcat.core.structures import (
    BioBank,
    CohortCfg,
    DiscoveryConfig,
    MetaAnalysisConfig,
    RunEnv,
    TraitType,
    Transform,
    make_conv,
)
from fedcat.platforms.base import JobSpec, PathSyntaxError, Platform
from fedcat.platforms.gcs_cromwell import gcs

# ── Config builders ───────────────────────────────────────────────────────────

AOU_WDLS = {
    "regenie_step1": "regenie_step1_scatter_aou",
    "regenie_step2_singlevar": "regenie_step2_single_variants_aou",
    "regenie_step2_aggregate": "regenie_step2_aggregates_aou",
    "conditioning_singlevar": "regenie_conditioning_singles_aou",
    "conditioning_aggregate": "regenie_conditioning_agg_aou",
}
_COHORT_SIZES = {"aou": {"eur": 1000, "afr": 800}, "ukb": {"eur": 2000}}


def cohort_raw(biobank: str, strata=("eur",), **overrides) -> dict:
    """A cohort as written in a config file: AoU with the repo's WDLs, UKB with fake applet and file ids."""
    if biobank == "aou":
        raw = {
            "workflows": {k: f"aou_regenie_workflows/wdls/{v}.wdl" for k, v in AOU_WDLS.items()},
            "pheno_file": "bucket/path/aou_pheno.tsv.gz",
            "cov_file": "bucket/path/aou_cov.tsv.gz",
            "cond_file": "bucket/path/aou_cond.tsv",
            "cat_vars": "seq_center",
            "pfile_template": "aou_pfiles/chr{chr}_v8",
        }
    else:
        raw = {
            "workflows": {k: f"applet-FAKE{i}" for i, k in enumerate(AOU_WDLS, 1)},
            "pheno_file": "file-AAAA1111",
            "cov_file": "file-BBBB2222",
            "cond_file": "file-CCCC3333",
            "cat_vars": "sequencing_centre",
            "pfile_template": "ukb_pfiles/chr{chr}_v3",
        }
    raw |= {
        "strata": list(strata),
        "step1_directory": f"{biobank}_step1",
        "root_results_directory": f"{biobank}_results",
        "step1_input_template": f"{biobank}_step1_input/{{stratum}}_step1",
        "regenie_annotations": [
            {
                "name": "annot1",
                **{
                    f"{kind}_file": f"{biobank}_annotations/{kind}s_chr{{chr}}.tsv.zst"
                    for kind in ("mask", "annot", "setlist")
                },
            }
        ],
        "cohort_sizes": {s: n for s, n in _COHORT_SIZES[biobank].items() if s in strata},
    }
    return raw | overrides


def make_config(cohorts: dict[str, dict] | None = None, **overrides) -> DiscoveryConfig:
    """A DiscoveryConfig for *cohorts* (raw, as from :func:`cohort_raw`; default AoU eur+afr and UKB eur)."""
    if cohorts is None:
        cohorts = {"aou": cohort_raw("aou", ["eur", "afr"]), "ukb": cohort_raw("ukb")}
    chroms = overrides.pop("chromosomes_to_run", ["1"])
    fields = {
        "meta_analysis_config": MetaAnalysisConfig(
            location=RunEnv.AOU,
            meta_workflow="aou_regenie_workflows/wdls/meta_analysis_discovery_aou.wdl",
            subset_workflow="aou_regenie_workflows/wdls/step0_subsetting_only_aou.wdl",
            schema_directory=overrides.pop("schema_directory", "wgs_annotations/v2026.2"),
        ),
        "local_env": RunEnv.AOU,
        "cohorts": make_conv().structure(cohorts, dict[BioBank, CohortCfg]),
        "pheno_names": ["TEST_PHENO"],
        "trait_type": TraitType.CONTINUOUS,
        "trait_transform": Transform.RINT,
        "chromosomes_to_run": chroms,
        "chromosomes_no_pgs": list(chroms),
        "local_run_path": "/tmp/test_cond_base",
        "run_prefix": "testrun",
        "cond_suffix": "testsuffix",
        "pgs_prefix": "XX_SCORE_CHR",
        "aou_bucket": "gs://dummy-bucket",
    }
    return DiscoveryConfig(**(fields | overrides))  # type: ignore[arg-type]


@pytest.fixture(autouse=True)
def _no_real_gcs(monkeypatch):
    """Keep tests off real storage (and Google's credential lookup, which waits to time out).

    A test that needs a client patches ``gcs.storage.Client`` with its own fake.
    """

    def _refuse(*_args, **_kwargs):
        raise RuntimeError("tests do not reach real GCS")

    monkeypatch.setattr(gcs.storage, "Client", _refuse)


@pytest.fixture(autouse=True)
def _fresh_async_primitives(monkeypatch):
    """Give each test its own module-level locks and semaphores.

    One waited on is bound to that event loop, and each test's ``asyncio.run`` makes a
    new one; fedcat itself runs a single loop per process.
    """
    import asyncio

    from fedcat.core import logging as fedcat_logging
    from fedcat.platforms.dnanexus import dx
    from fedcat.platforms.dnanexus import platform as dx_platform
    from fedcat.platforms.gcs_cromwell import cromwell

    for module in (fedcat_logging, dx, dx_platform, cromwell):
        for name, value in list(vars(module).items()):
            if isinstance(value, (asyncio.Lock, asyncio.Semaphore)):
                # sized as at import, so a permit one test leaks is not carried into the next
                size = _PRIMITIVE_SIZES.setdefault((module.__name__, name), getattr(value, "_value", 0))
                fresh = type(value)() if isinstance(value, asyncio.Lock) else type(value)(size)
                monkeypatch.setattr(module, name, fresh)


_PRIMITIVE_SIZES: dict[tuple[str, str], int] = {}


class FakePlatform(Platform):
    """In-memory :class:`Platform` for exercising engine/job/conditioning code
    without touching GCS or DNAnexus.

    ``existing`` is the set of refs :meth:`exists` reports as present. Submitted
    jobs are recorded in ``submitted`` and always "succeed" when waited on.
    """

    def __init__(
        self,
        kind: str = "fake",
        existing: set[str] | None = None,
        missing_targets: set[str] | None = None,
        connect_error: str | None = None,
        submit_connect_error: str | None = None,
        exists_fn=None,
    ) -> None:
        self.kind = kind
        self.existing: set[str] = set(existing or ())
        self.missing_targets: set[str] = set(missing_targets or ())
        self.connect_error = connect_error
        self.submit_connect_error = submit_connect_error  # raised only when need_submit=True
        self.exists_fn = exists_fn  # optional Callable[[str], bool] overriding `existing`
        self.exists_seen: list[str] = []
        self.submitted: list[JobSpec] = []
        self.uploaded: list[tuple[Path, str]] = []
        self.connected = False
        self.aborted: list[str] = []
        self.waited: list[tuple[str, float | None]] = []  # (job id, initial_delay)
        self.abort_error: str | None = None  # when set, abort() raises AbortError with this message
        self.task_statuses: dict[str, list] = {}  # job id -> [TaskStatus, ...]
        self.reported_cost = None  # Cost returned by job_costs()
        self.hourly_cost: float | None = None  # returned by task_hourly_cost()
        self.files: dict[str, dict[str, str]] = {}  # folder -> {file name: text} served by download()
        self.listing: dict[str, list[str]] = {}  # folder -> file names served by list_files()

    async def connect(self, *, need_submit: bool) -> None:
        from fedcat.platforms.base import AuthError

        if self.connect_error:
            raise AuthError(self.connect_error)
        if need_submit and self.submit_connect_error:
            raise AuthError(self.submit_connect_error)
        self.connected = True

    def canonical_data_ref(self, ref: str) -> str:
        return ref

    def canonical_path_ref(self, ref: str) -> str:
        return ref

    async def exists(self, ref: str) -> bool:
        if not ref:
            raise PathSyntaxError("empty reference")
        self.exists_seen.append(ref)
        if self.exists_fn is not None:
            return self.exists_fn(ref)
        return ref in self.existing

    async def list_files(self, folder: str, pattern: str | None = None) -> list[str]:
        return [n for n in self.listing.get(folder, []) if pattern is None or re.match(pattern, n)]

    async def download(self, folder, pattern, dest: Path, max_workers=None) -> None:
        dest.mkdir(parents=True, exist_ok=True)
        for name, text in self.files.get(folder, {}).items():
            if re.match(pattern, name):
                (dest / name).write_text(text)

    async def upload(self, file: Path, destination: str, *, overwrite: bool = False) -> str:
        self.uploaded.append((file, destination))
        return destination

    async def check_targets(self, targets: list[str]) -> list[str]:
        return [t for t in targets if t in self.missing_targets]

    async def submit(self, spec: JobSpec) -> str:
        self.submitted.append(spec)
        return f"fake-job-{len(self.submitted)}"

    async def wait(
        self,
        job_id: str,
        *,
        poll_interval: float,
        description: str,
        initial_delay: float | None = None,
    ) -> tuple[str, bool]:
        self.waited.append((job_id, initial_delay))
        return "Succeeded", True

    async def abort(self, job_id: str) -> None:
        if self.abort_error:
            from fedcat.platforms.base import AbortError

            raise AbortError(self.abort_error)
        self.aborted.append(job_id)

    async def task_counts(self, job_ids: list[str]) -> Counter:
        if self.submit_connect_error:  # no job service, as when a Cromwell server is down
            from fedcat.platforms.base import TrackingError

            raise TrackingError(self.submit_connect_error)
        return Counter(s for job_id in job_ids for s in self.task_statuses.get(job_id, []))

    async def job_costs(self, job_ids: list[str]):
        return self.reported_cost

    def task_hourly_cost(self, target: str) -> float | None:
        return self.hourly_cost

    def results_root(self, root_results_directory: str) -> str:
        return f"/{root_results_directory.lstrip('/')}"


@pytest.fixture
def fake_platform() -> FakePlatform:
    """A fresh FakePlatform."""
    return FakePlatform()
