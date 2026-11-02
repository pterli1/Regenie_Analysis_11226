"""Integration tests running parts of the pipeline together: phases, conditioning and staging."""

from __future__ import annotations

import asyncio
import io
import re
import shutil
import tarfile
import weakref
from pathlib import Path

import polars as pl
import pytest

from fedcat.core.conditioning import ConditioningStep
from fedcat.core.engine import _merge_step0_subsets
from fedcat.core.structures import BioBank, CondType, DiscoveryConfig
from fedcat.stats.subset_results import _PARQUET_SCHEMA_FILES, subset_results_for_phenotype
from fedcat.tests.conftest import FakePlatform, cohort_raw, make_config


def _make_aou_cfg(tmp_path: Path, *, schema_directory: str = "wgs_schema/chr{chr}") -> DiscoveryConfig:
    """AoU eur only, on chromosomes 1 and X, with one file as both phenotype and covariate file."""
    aou = cohort_raw("aou", pheno_file="bucket/shared_input.tsv.gz", cov_file="bucket/shared_input.tsv.gz")
    return make_config(
        {"aou": aou},
        chromosomes_to_run=["1", "X"],
        local_run_path=str(tmp_path / "cond_base"),
        schema_directory=schema_directory,
    )


def _remote_findings(cfg, platform, phases_to_run=None):
    """Run the remote-files check for a single-AoU-platform config."""
    from fedcat.core.checks import CheckContext, check_remote_files
    from fedcat.core.structures import Phase

    ctx = CheckContext(
        config=cfg,
        platforms={BioBank.AOU: platform},
        # Annotation/schema references are only checked for the phases that use
        # them, so ask for the aggregate + meta phases by default.
        phases_to_run=[Phase.STEP2_AGG, Phase.META] if phases_to_run is None else phases_to_run,
    )
    return asyncio.run(check_remote_files(ctx))


def test_validate_remote_files_checks_meta_schema(tmp_path):

    cfg = _make_aou_cfg(tmp_path)
    platform = FakePlatform(kind="gcs_cromwell", exists_fn=lambda _p: True)

    findings = _remote_findings(cfg, platform)
    assert findings == []

    # The meta schema files are among the references that got checked.
    for chrom in cfg.chromosomes_to_run:
        for schema_file in _PARQUET_SCHEMA_FILES:
            assert any(seen.endswith(f"chr{chrom}/{schema_file}") for seen in platform.exists_seen)


@pytest.mark.parametrize(
    ("present_exts", "expect_ok"),
    [
        ((".bed", ".bim", ".fam"), True),
        ((".pgen", ".pvar", ".psam"), True),
        ((".bed", ".bim", ".fam", ".pgen", ".pvar", ".psam"), True),
        ((".bed", ".bim"), False),  # incomplete PLINK1 set
        ((".bed", ".bim", ".pgen", ".pvar"), False),  # neither set complete
        ((), False),
    ],
)
def test_validate_remote_files_step1_accepts_either_fileset(tmp_path, present_exts, expect_ok):

    cfg = _make_aou_cfg(tmp_path)
    step1_prefix = "aou_step1_input/eur_step1"

    def _exists(path: str) -> bool:
        p = path.removeprefix("gs://")
        if p.startswith(step1_prefix):
            return p[len(step1_prefix) :] in present_exts
        return True

    findings = _remote_findings(cfg, FakePlatform(kind="gcs_cromwell", exists_fn=_exists))
    step1_findings = [f for f in findings if "step1_input" in f.label]
    assert (step1_findings == []) is expect_ok


def _make_disc_run(cfg):
    """Return a discovery run for `cfg`."""
    from fedcat.core.engine import DiscoveryRun
    from fedcat.core.structures import Phase

    run_dir = Path(cfg.local_run_path) / f"{cfg.run_prefix}_discovery_run"
    return DiscoveryRun(
        config=cfg,
        run_path=run_dir,
        variant_density_blockfile=Path("variant_density_blocks.tsv"),
        phases_to_run=[Phase.COND_SINGLES],
        current_phase=Phase.COND_SINGLES,
    )


def test_cond_tarball_override_rejected_for_multi_pheno(tmp_path):
    """A direct step-0 tarball names one pheno's subset; refuse a multi-pheno run."""
    cfg = _make_aou_cfg(tmp_path)
    cfg.pheno_names = ["PHENO_A", "PHENO_B"]
    disc = _make_disc_run(cfg)
    disc.cond_singles_tarball = tmp_path / "external_singles.tar.gz"

    from fedcat.core.engine import CondOutcome

    assert asyncio.run(disc._run_conditioning(CondType.SINGLEVAR)) is CondOutcome.FAILED


def test_cond_tarball_override_flows_into_fresh_step(monkeypatch, tmp_path):
    """With a single pheno, the override tarball becomes the fresh step's input path."""
    cfg = _make_aou_cfg(tmp_path)
    disc = _make_disc_run(cfg)
    override = tmp_path / "external_singles.tar.gz"
    override.write_bytes(b"")
    disc.cond_singles_tarball = override

    captured: dict[str, Path | None] = {}

    async def _fake_run_to_completion(self):
        captured["tarball"] = self.input_tarball_path
        # Signal "no more signals" so the outer loop stops after this step.
        self.required = False
        return True

    # Stop the state machine right after construction so we only inspect wiring.
    monkeypatch.setattr(ConditioningStep, "run_to_completion", _fake_run_to_completion)
    monkeypatch.setattr(ConditioningStep, "persist_state", _noop_persist)

    async def _drive():
        from fedcat.core.logging import STATE_WRITE_QUEUE

        STATE_WRITE_QUEUE.start()
        try:
            await disc._run_conditioning(CondType.SINGLEVAR)
        finally:
            await STATE_WRITE_QUEUE.stop()

    asyncio.run(_drive())
    assert captured["tarball"] == override


def test_conditioning_pause_is_passed_on_not_recorded_as_failure(monkeypatch, tmp_path):
    """Ctrl-C mid-conditioning pauses: the cancellation propagates and no phenotype is marked failed."""
    cfg = _make_aou_cfg(tmp_path)
    disc = _make_disc_run(cfg)
    started = asyncio.Event()

    async def _run_forever(self):
        started.set()
        await asyncio.Event().wait()

    monkeypatch.setattr(ConditioningStep, "run_to_completion", _run_forever)
    monkeypatch.setattr(ConditioningStep, "persist_state", _noop_persist)

    async def _drive():
        task = asyncio.create_task(disc._run_conditioning(CondType.SINGLEVAR))
        await started.wait()
        task.cancel()
        with pytest.raises(asyncio.CancelledError):
            await task

    asyncio.run(_drive())
    assert disc.cond_status.get(CondType.SINGLEVAR, {}) == {}


@pytest.mark.parametrize(
    "step_state, job_states, expected_step, expected_jobs",
    [
        # Paused mid-run with nothing failed: resume exactly as it was.
        ("RUN_JOBS", ["RUNNING", "COMPLETED"], "RUN_JOBS", ["RUNNING", "COMPLETED"]),
        ("SORT_INPUT", [], "SORT_INPUT", []),
        # Failed jobs are resubmitted; the others are left alone.
        ("ERROR", ["ERROR", "COMPLETED"], "RUN_JOBS", ["EMPTY", "COMPLETED"]),
        ("RUN_JOBS", ["ERROR", "RUNNING"], "RUN_JOBS", ["EMPTY", "RUNNING"]),
    ],
)
def test_force_job_retry_only_resets_what_failed(tmp_path, step_state, job_states, expected_step, expected_jobs):
    from fedcat.core.engine import DiscoveryRun
    from fedcat.core.job import Job
    from fedcat.core.structures import JobState, StepState

    step = ConditioningStep.from_discovery_config(_make_aou_cfg(tmp_path), "TEST_PHENO", CondType.SINGLEVAR)
    step.state = StepState[step_state]
    step.jobs = [Job(biobank=BioBank.AOU, description=f"job{i}") for i in range(len(job_states))]
    for job, state in zip(step.jobs, job_states):
        job.state = JobState[state]

    DiscoveryRun._retry_failed_conditioning(step)

    assert step.state is StepState[expected_step]
    assert [j.state for j in step.jobs] == [JobState[s] for s in expected_jobs]


def test_force_job_retry_reruns_stage_a_step_failed_in(tmp_path):
    from fedcat.core.engine import DiscoveryRun
    from fedcat.core.structures import StepState

    step = ConditioningStep.from_discovery_config(_make_aou_cfg(tmp_path), "TEST_PHENO", CondType.SINGLEVAR)
    step.state = StepState.CHECK_RESULTS
    step.set_state(StepState.ERROR)
    DiscoveryRun._retry_failed_conditioning(step)
    assert step.state is StepState.CHECK_RESULTS


def _noop_persist(*_args, **_kwargs):
    """Do nothing, in place of saving state."""


def test_merge_step0_subsets_raises_for_missing_flat_annotations(tmp_path):
    pheno = "TEST_PHENO"
    pheno_dl = tmp_path / "parts"
    pheno_dl.mkdir()
    source_payload = tmp_path / "payload"
    source_payload.mkdir()
    (source_payload / "tophits_input").write_text("ID\t6\tNA\tNA\tNA\tBURDEN\t1\t100\t200\n")

    with tarfile.open(pheno_dl / f"{pheno}_agg_cond_step0_subset_chr1.tar.gz", "w:gz") as tf:
        tf.add(source_payload / "tophits_input", arcname="tophits_input")

    with pytest.raises(ValueError, match="Missing flat annotation files"):
        _merge_step0_subsets(pheno, pheno_dl, tmp_path)


def _singles_tarball(path: Path, variant_id: str, chrom: str) -> Path:
    """Write a singles step-0 subset tarball of one variant."""
    payload = path.parent / f"{path.name}.tophits"
    payload.write_text(f"{variant_id}\t7\t0.1\t0.2\t0.01\tADD\t{chrom}\t10\t10\n")
    with tarfile.open(path, "w:gz") as tf:
        tf.add(payload, arcname="tophits_input")
    return path


def _input_names(spec) -> dict:
    """A submitted job's inputs keyed by input name, without the workflow prefix."""
    return {k.rsplit(".", 1)[-1]: v for k, v in spec.inputs.items()}


class _Step0Platform(FakePlatform):
    """Serves tarballs by folder, and marks each submitted subset task complete."""

    def __init__(self) -> None:
        super().__init__(kind="gcs_cromwell")
        self.tarballs: dict[str, list[Path]] = {}

    async def submit(self, spec):
        job_id = await super().submit(spec)
        inputs = _input_names(spec)
        folder = inputs["output_folder"].removeprefix("gs://")
        self.listing.setdefault(folder, []).extend(
            f"testrun_step0_subset_chr{c}.log" for c in inputs["chromosomes"]
        )
        return job_id

    async def download(self, folder, pattern, dest: Path, max_workers=None) -> None:
        dest.mkdir(parents=True, exist_ok=True)
        for src in self.tarballs.get(folder, []):
            if re.match(pattern, src.name):
                shutil.copy(src, dest / src.name)


def _drive_phase(disc) -> None:
    """Run the current phase; callers patch out ``persist_state``, so no state write queue is needed."""
    asyncio.run(disc._run_phase(reset_jobs=None, reset_jobs_list=[]))


def _make_step0_run(tmp_path, monkeypatch):
    """Return an AoU run at the step0_subset phase, with two phenotypes and a fake platform."""
    from fedcat.core.structures import Phase

    cfg = _make_aou_cfg(tmp_path)
    cfg.pheno_names = ["TEST_PHENO", "TEST_PHENO_F"]
    disc = _make_disc_run(cfg)
    disc.phases_to_run = [Phase.STEP0_SUBSET]
    disc.current_phase = Phase.STEP0_SUBSET
    platform = _Step0Platform()
    disc.platforms = {BioBank.AOU: platform}
    monkeypatch.setattr(type(disc), "persist_state", _noop_persist)
    return disc, platform


def test_step0_subset_runs_missing_chromosomes_then_merges(monkeypatch, tmp_path):
    disc, platform = _make_step0_run(tmp_path, monkeypatch)
    singles_folder = disc._step0_folder(CondType.SINGLEVAR)
    agg_folder = disc._step0_folder(CondType.AGGREGATE)
    assert singles_folder.endswith("all_biobank_meta/step0_subsets/singles_p1e-06")
    platform.listing[singles_folder] = ["testrun_step0_subset_chr1.log"]
    src = tmp_path / "src"
    src.mkdir()
    platform.tarballs[singles_folder] = [
        _singles_tarball(src / "TEST_PHENO_singles_cond_step0_subset_chr1.tar.gz", "v1", "1"),
        _singles_tarball(src / "TEST_PHENO_singles_cond_step0_subset_chrX.tar.gz", "vX", "X"),
        # A phenotype whose name extends TEST_PHENO must not be merged into it.
        _singles_tarball(src / "TEST_PHENO_F_singles_cond_step0_subset_chr1.tar.gz", "vF", "1"),
    ]
    stale_agg = disc.run_path / "TEST_PHENO_agg_cond_step0_subset.tar.gz"
    stale_agg.write_bytes(b"old")

    _drive_phase(disc)

    assert not disc.had_error
    jobs = {i["run_uncond_singles"]: i for i in map(_input_names, platform.submitted)}
    assert jobs[True]["chromosomes"] == ["X"]
    assert jobs[False]["chromosomes"] == ["1", "X"]
    assert jobs[False]["output_folder"] == f"gs://{agg_folder}"
    with tarfile.open(disc.run_path / "TEST_PHENO_singles_cond_step0_subset.tar.gz") as tf:
        ids = sorted(line.split("\t")[0] for line in tf.extractfile("tophits_input").read().decode().splitlines())
    assert ids == ["v1", "vX"]
    assert not stale_agg.exists()

    # New thresholds subset afresh into their own folder; the unchanged ones are left alone.
    disc.config.cond_p_singles_retest = 5e-8
    _drive_phase(disc)
    assert not disc.had_error
    new_jobs = [_input_names(spec) for spec in platform.submitted[2:]]
    assert [i["output_folder"] for i in new_jobs] == [f"gs://{disc._step0_folder(CondType.SINGLEVAR)}"]
    assert disc._step0_folder(CondType.SINGLEVAR).endswith("singles_p5e-08")
    assert new_jobs[0]["chromosomes"] == ["1", "X"]


def test_meta_phase_resumes_jobs_without_restaging(monkeypatch, tmp_path):
    from fedcat.core.job import Job
    from fedcat.core.structures import JobState, Phase

    cfg = _make_aou_cfg(tmp_path)
    disc = _make_disc_run(cfg)
    disc.phases_to_run = [Phase.META]
    disc.current_phase = Phase.META
    disc.platforms = {BioBank.AOU: FakePlatform(kind="gcs_cromwell")}
    monkeypatch.setattr(type(disc), "persist_state", _noop_persist)
    staged: list[bool] = []

    async def _stage(_self):
        staged.append(True)
        return True

    monkeypatch.setattr(type(disc), "_stage_cohort_results_to_meta", _stage)

    running = Job(biobank=BioBank.AOU, description="meta", recipe={"op": "meta"}, id="wf-1")
    running.state = JobState.RUNNING
    running._parent = weakref.ref(disc)
    disc.jobs = [running]
    _drive_phase(disc)
    assert staged == []
    assert running.state == JobState.COMPLETED

    disc.jobs = []
    _drive_phase(disc)
    assert staged == [True]
    assert [j.state for j in disc.jobs] == [JobState.COMPLETED]


def test_meta_analyse_aggs_preserves_add_chr_burden_ids(monkeypatch, tmp_path):
    """Regression guard for the ADD_chr parsing bug in aggregate conditioning merges."""
    cfg = _make_aou_cfg(tmp_path)
    # Keep a single chromosome for a compact integration fixture.
    cfg.chromosomes_to_run = ["1"]

    step = ConditioningStep.from_discovery_config(cfg, "TEST_PHENO", CondType.AGGREGATE)
    step_path = step.step_path
    step_path.mkdir(parents=True, exist_ok=True)
    (step_path / "regenie_output_aou").mkdir(parents=True, exist_ok=True)

    # Make glob discovery include all aggregate result types.
    for fname in [
        "chunk_eur_chr1_x_burden_TEST_PHENO.regenie.gz",
        "chunk_eur_chr1_x_mac10_TEST_PHENO.regenie.gz",
        "chunk_eur_chr1_singlevars_TEST_PHENO.regenie.gz",
    ]:
        (step_path / "regenie_output_aou" / fname).write_text("dummy")

    # Required annotation inputs for ACAT branch selection.
    step0_path = step.run_path / "agg_cond_step0"
    step0_path.mkdir(parents=True, exist_ok=True)
    (step0_path / "annotations").write_text("")
    (step0_path / "masks").write_text("")

    hit_id = "GENE1.MASK1_ADD_chr1"
    (step_path / "tophits_input").write_text(f"{hit_id}\t8\tNA\tNA\tNA\tBURDEN\t1\t100\t200\n")

    import fedcat.core.conditioning as cond_mod

    async def _fake_meta(_cohorts, output_file: Path, res_type: str, **_kwargs):
        output_file.parent.mkdir(parents=True, exist_ok=True)
        header = "CHROM\tPOS\tID\tA0\tA1\tA1FREQ\tN\tTEST\tBETA\tSE\tZ\tLOG10P\tQ\tQ_P\tI2\tN_STUDIES"
        if res_type == "burden":
            row = f"1\t100\t{hit_id}\tA\tG\t0.01\t1000\tBURDEN\t0.2\t0.05\t4\t12\t0\t1\t0\t1"
        elif res_type == "mac10":
            row = "1\t110\tGENE2.MASK2_ADD_chr1\tA\tG\t0.01\t1000\tBURDEN\t0.1\t0.05\t2\t4\t0\t1\t0\t1"
        else:
            row = "1\t120\t1:120:A:G\tA\tG\t0.01\t1000\tADD\t0.1\t0.05\t2\t4\t0\t1\t0\t1"
        output_file.write_text(f"{header}\n{row}\n")

    def _fake_compute_acat(**kwargs):
        out = Path(kwargs["output_file"])
        out.write_text(
            "CHR\tSTART\tEND\tID\tUNIT\tMASK\tTHRESHOLD\tN_AGG_COMBINED\tN_SINGLEVAR_COMBINED\tLOG10P_ACAT\n"
        )

    async def _noop_async(*_args, **_kwargs):
        return None

    original_pl_read_csv = cond_mod.pl.read_csv

    def _patched_pl_read_csv(source, **kwargs):
        # Test files use .zst extension but are plain text; bypass extension-based decompression.
        if isinstance(source, (str, Path)) and str(source).endswith(".zst"):
            source = io.BytesIO(Path(source).read_bytes())
        return original_pl_read_csv(source, **kwargs)

    monkeypatch.setattr(cond_mod, "meta_analyse_regenie_async", _fake_meta)
    monkeypatch.setattr(cond_mod, "compute_acat", _fake_compute_acat)
    monkeypatch.setattr(cond_mod.pl, "read_csv", _patched_pl_read_csv)
    monkeypatch.setattr(ConditioningStep, "_upload_meta_outputs_to_cloud", _noop_async)
    monkeypatch.setattr(ConditioningStep, "persist_state", _noop_persist)

    asyncio.run(step._meta_analyse_aggs())

    summary = pl.read_csv(step_path / "summary_output", separator="\t", null_values="NA")
    assert summary.height == 1
    assert summary[0, "ID"] == hit_id
    assert summary[0, "TEST"] == "BURDEN"
    assert summary[0, "LOG10P_OUT"] == pytest.approx(12.0)


def _run_meta_analyse_aggs(monkeypatch, tmp_path, result_files: list[str], acat_singles_stratum_rare: bool):
    """Drive ``_meta_analyse_aggs`` on fake per-stratum outputs; returns the meta and ACAT calls."""
    import fedcat.core.conditioning as cond_mod

    cfg = _make_aou_cfg(tmp_path)
    cfg.chromosomes_to_run = ["1"]
    cfg.aggregate_freqs = ["0.001", "0.01"]
    cfg.acat_singles_stratum_rare = acat_singles_stratum_rare
    step = ConditioningStep.from_discovery_config(cfg, "TEST_PHENO", CondType.AGGREGATE)
    results = step.step_path / "regenie_output_aou"
    results.mkdir(parents=True, exist_ok=True)
    for fname in result_files:
        (results / fname).write_text("dummy")
    step0_path = step.run_path / "agg_cond_step0"
    step0_path.mkdir(parents=True, exist_ok=True)
    (step0_path / "annotations").write_text("")
    (step0_path / "masks").write_text("")
    (step.step_path / "tophits_input").write_text("GENE1.MASK1.0.001\t8\tNA\tNA\tNA\tACAT\t1\t100\t200\n")

    meta_calls: list[tuple[str, str, object]] = []
    acat_calls: list[dict] = []
    header = "CHROM\tPOS\tID\tA0\tA1\tA1FREQ\tN\tTEST\tBETA\tSE\tZ\tLOG10P\tQ\tQ_P\tI2\tN_STUDIES\n"

    async def _fake_meta(_cohorts, output_file: Path, res_type: str, **kwargs):
        meta_calls.append((output_file.name, res_type, kwargs.get("aaf_bins")))
        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_text(
            header + "1\t120\t1:120:A:G\tA\tG\t0.0005\t1000\tADD\t0.1\t0.05\t2\t4\t0\t1\t0\t1\n"
        )

    def _fake_compute_acat(**kwargs):
        acat_calls.append(kwargs)
        assert Path(kwargs["single_variants_file"]).exists()
        Path(kwargs["output_file"]).write_text(
            "CHR\tSTART\tEND\tID\tUNIT\tMASK\tTHRESHOLD\tN_AGG_COMBINED\tN_SINGLEVAR_COMBINED\tLOG10P_ACAT\n"
            "1\t100\t200\tGENE1.MASK1.0.001\tGENE1\tMASK1\t0.001\t0\t1\t6\n"
        )

    async def _noop_async(*_args, **_kwargs):
        return None

    original_pl_read_csv = cond_mod.pl.read_csv

    def _patched_pl_read_csv(source, **kwargs):
        if isinstance(source, (str, Path)) and str(source).endswith(".zst"):
            source = io.BytesIO(Path(source).read_bytes())
        return original_pl_read_csv(source, **kwargs)

    monkeypatch.setattr(cond_mod, "meta_analyse_regenie_async", _fake_meta)
    monkeypatch.setattr(cond_mod, "compute_acat", _fake_compute_acat)
    monkeypatch.setattr(cond_mod.pl, "read_csv", _patched_pl_read_csv)
    monkeypatch.setattr(ConditioningStep, "_upload_meta_outputs_to_cloud", _noop_async)
    monkeypatch.setattr(ConditioningStep, "persist_state", _noop_persist)
    asyncio.run(step._meta_analyse_aggs())
    return step, meta_calls, acat_calls


_ALL_BUILDS = [
    "chunk_eur_chr1_x_burden_TEST_PHENO.regenie.gz",
    "chunk_eur_chr1_x_mac10_TEST_PHENO.regenie.gz",
    "chunk_eur_chr1_singlevars_TEST_PHENO.regenie.gz",
]


def test_meta_analyse_aggs_acat_uses_per_stratum_binned_singles(monkeypatch, tmp_path):
    _, meta_calls, acat_calls = _run_meta_analyse_aggs(monkeypatch, tmp_path, _ALL_BUILDS, True)
    assert ("META_single_variants_acat.tsv.zst", "single_variants", ["0.001", "0.01"]) in meta_calls
    # the plain singles meta is still built, unbinned
    assert ("META_single_variants.tsv.zst", "single_variants", None) in meta_calls
    (acat,) = acat_calls
    assert Path(acat["single_variants_file"]).name == "META_single_variants_acat.tsv.zst"


def test_meta_analyse_aggs_legacy_rule_uses_plain_singles(monkeypatch, tmp_path):
    _, meta_calls, acat_calls = _run_meta_analyse_aggs(monkeypatch, tmp_path, _ALL_BUILDS, False)
    assert all(bins is None for _, _, bins in meta_calls)
    (acat,) = acat_calls
    assert Path(acat["single_variants_file"]).name == "META_single_variants.tsv.zst"


def test_meta_analyse_aggs_singles_only_step(monkeypatch, tmp_path):
    """Every mask below --minMAC in every stratum: no burden or mac10 meta, singles still ACAT'd."""
    step, _, acat_calls = _run_meta_analyse_aggs(
        monkeypatch, tmp_path, ["chunk_eur_chr1_singlevars_TEST_PHENO.regenie.gz"], True
    )
    (acat,) = acat_calls
    assert acat["burden_file"] is None
    summary = pl.read_csv(step.step_path / "summary_output", separator="\t", null_values="NA")
    assert summary[0, "LOG10P_OUT"] == pytest.approx(6.0)


def test_meta_analyse_aggs_burden_only_step(monkeypatch, tmp_path):
    """No stratum returned singles or mac10: ACAT gets an empty singles input, not a missing file."""
    _, _, acat_calls = _run_meta_analyse_aggs(
        monkeypatch, tmp_path, ["chunk_eur_chr1_x_burden_TEST_PHENO.regenie.gz"], True
    )
    (acat,) = acat_calls
    assert Path(acat["single_variants_file"]).read_text() == "ID\tA1FREQ\tLOG10P\n"
    assert acat["burden_file"] is None


def test_subset_parquet_schema_gives_the_kept_units_the_same_acat(tmp_path):
    """Step-0 subsetting with schema_directory cuts the meta schema to the selected units; their ACAT must not change."""
    from fedcat.stats.acat import compute_acat
    from fedcat.stats.subset_results import _subset_annotations_parquet
    from fedcat_rs import convert_tsv_to_parquet_rs

    annots = tmp_path / "annots.tsv"
    annots.write_text("1:100:A:G\tG\tpLoF\n1:200:A:G\tG\tmissense\n1:250:A:G\tH\tpLoF\n1:300:A:G\tH\tmissense\n")
    masks = tmp_path / "masks.tsv"
    masks.write_text("lof\tpLoF\nall\tpLoF,missense\n")
    full, subset = tmp_path / "full", tmp_path / "subset"
    convert_tsv_to_parquet_rs([str(annots)], [str(masks)], str(full), None)
    _subset_annotations_parquet(full, {"G"}, subset)
    assert sorted(p.name for p in subset.iterdir()) == sorted(_PARQUET_SCHEMA_FILES)

    singles = tmp_path / "singles.tsv"
    pl.DataFrame(
        {
            "ID": ["1:100:A:G", "1:200:A:G", "1:250:A:G", "1:300:A:G"],
            "A1FREQ": [0.0005] * 4,
            "LOG10P": [3.0, 5.0, 2.0, 4.0],
        }
    ).write_csv(singles, separator="\t")

    def _acat(schema: Path) -> pl.DataFrame:
        out = tmp_path / f"acat_{schema.name}.tsv"
        compute_acat(single_variants_file=singles, schema_dir=schema, output_file=out, n_threads=1)
        return pl.read_csv(out, separator="\t").sort("ID")

    expected = _acat(full).filter(pl.col("UNIT") == "G")
    assert expected.height == 2
    assert _acat(subset).equals(expected)


def test_subset_results_applies_separate_agg_and_singles_thresholds(monkeypatch, tmp_path):
    """Step-0 subsets must use cond_p_agg_retest and cond_p_singles_retest independently."""
    monkeypatch.chdir(tmp_path)
    acat = tmp_path / "acat.tsv"
    pl.DataFrame(
        {
            "ID": ["G1.M1", "G2.M1"],
            "UNIT": ["G1", "G2"],
            "CHR": ["1", "1"],
            "START": [100, 200],
            "END": [150, 250],
            "LOG10P_ACAT": [7.5, 6.5],
        }
    ).write_csv(acat, separator="\t")
    singles = tmp_path / "singles.tsv"
    pl.DataFrame({"CHROM": ["1", "1"], "POS": [10, 20], "ID": ["v1", "v2"], "LOG10P": [7.5, 6.5]}).write_csv(
        singles, separator="\t"
    )
    schema_dir = tmp_path / "schema"
    schema_dir.mkdir()
    for name in _PARQUET_SCHEMA_FILES:
        (schema_dir / name).touch()
    monkeypatch.setattr("fedcat.stats.subset_results._subset_annotations_parquet", lambda *a, **k: None)

    subset_results_for_phenotype(
        "P", acat_path=acat, singles_path=singles, schema_dir=schema_dir, agg_threshold=7.0, singles_threshold=6.0
    )

    def _tophit_ids(tar_name: str) -> list[str]:
        with tarfile.open(tmp_path / tar_name) as tar:
            member = tar.extractfile("tophits_input")
            assert member is not None
            return [line.split("\t")[0] for line in member.read().decode().splitlines()]

    assert _tophit_ids("P_agg_cond_step0_subset.tar.gz") == ["G1.M1"]
    assert _tophit_ids("P_singles_cond_step0_subset.tar.gz") == ["v1", "v2"]


def test_subset_results_keeps_both_tests_and_dotted_units(monkeypatch, tmp_path):
    """ACAT and burden share UNIT.MASK.BIN IDs: both rows are kept, and units may contain dots."""
    monkeypatch.chdir(tmp_path)
    acat = tmp_path / "acat.tsv"
    pl.DataFrame(
        {
            "ID": ["ENST1.2.M1.0.001"],
            "UNIT": ["ENST1.2"],
            "CHR": ["1"],
            "START": [100],
            "END": [150],
            "LOG10P_ACAT": [8.0],
        }
    ).write_csv(acat, separator="\t")
    burden = tmp_path / "burden.tsv"
    pl.DataFrame(
        {
            "CHROM": ["1"],
            "POS": [100],
            "ID": ["ENST1.2.M1.0.001"],
            "A1FREQ": [0.001],
            "BETA": [0.5],
            "SE": [0.1],
            "LOG10P": [7.5],
        }
    ).write_csv(burden, separator="\t")
    schema_dir = tmp_path / "schema"
    schema_dir.mkdir()
    for name in _PARQUET_SCHEMA_FILES:
        (schema_dir / name).touch()
    subset_units: list[set[str]] = []
    monkeypatch.setattr(
        "fedcat.stats.subset_results._subset_annotations_parquet",
        lambda _schema, gene_set, _out: subset_units.append(gene_set),
    )

    subset_results_for_phenotype(
        "P", acat_path=acat, burden_path=burden, schema_dir=schema_dir, agg_threshold=7.0
    )

    with tarfile.open(tmp_path / "P_agg_cond_step0_subset.tar.gz") as tar:
        member = tar.extractfile("tophits_input")
        assert member is not None
        rows = [line.split("\t") for line in member.read().decode().splitlines()]
    assert [(r[0], r[5]) for r in rows] == [("ENST1.2.M1.0.001", "ACAT"), ("ENST1.2.M1.0.001", "BURDEN")]
    # burden row gets the unit's ACAT bounds
    assert rows[1][7:9] == ["100", "150"]
    assert subset_units == [{"ENST1.2"}]


@pytest.mark.parametrize("one_step", [False, True])
def test_conditioning_stopped_early_is_incomplete_not_failed(monkeypatch, tmp_path, one_step):
    """cond_max_steps writes partial results and --cond-one-step stops: neither is an error."""
    from fedcat.core.engine import CondOutcome

    cfg = _make_aou_cfg(tmp_path)
    cfg.cond_max_steps = 2
    disc = _make_disc_run(cfg)
    disc.cond_one_step = one_step
    written: list[bool] = []

    async def _step_done(self):
        self.step_number = 2
        return True

    async def _write(self, partial=False):
        written.append(partial)
        return True

    monkeypatch.setattr(ConditioningStep, "run_to_completion", _step_done)
    monkeypatch.setattr(ConditioningStep, "write_final_results", _write)
    monkeypatch.setattr(ConditioningStep, "persist_state", _noop_persist)

    async def _drive():
        from fedcat.core.logging import STATE_WRITE_QUEUE

        STATE_WRITE_QUEUE.start()
        try:
            await disc._run_phase(reset_jobs=None, reset_jobs_list=[], attempt_job_reset=False)
        finally:
            await STATE_WRITE_QUEUE.stop()

    asyncio.run(_drive())
    assert not disc.had_error
    assert disc._stop_requested is one_step
    assert written == ([] if one_step else [True])
    # incomplete phenotypes are resumed, not skipped, next time
    assert disc.cond_status[CondType.SINGLEVAR] == {cfg.pheno_names[0]: CondOutcome.INCOMPLETE}


def test_meta_analyse_singles_keeps_variant_ids_as_in_the_pvar(monkeypatch, tmp_path):
    """IDs pass through unchanged, so the next step's jobs can find them in the pvar; CHROM is normalised."""
    import fedcat.core.conditioning as cond_mod

    cfg = _make_aou_cfg(tmp_path)
    cfg.chromosomes_to_run = ["1"]
    cfg.cond_check_coll = False
    step = ConditioningStep.from_discovery_config(cfg, "TEST_PHENO", CondType.SINGLEVAR)
    results = step.step_path / "regenie_output_aou"
    results.mkdir(parents=True)
    (results / "x_eur_chr1_singlevar_TEST_PHENO.regenie.gz").write_text("dummy")
    (step.step_path / "tophits_input").write_text("chr1:120:A:G\t8\t0.01\t0.1\t0.02\tADD\tchr1\t120\t120\n")

    async def _fake_meta(_cohorts, output_file: Path, _res_type: str, **_kwargs):
        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_text(
            "CHROM\tPOS\tID\tA0\tA1\tA1FREQ\tN\tTEST\tBETA\tSE\tZ\tLOG10P\tQ\tQ_P\tI2\tN_STUDIES\n"
            "1\t120\tchr1:120:A:G\tA\tG\t0.01\t1000\tADD\t0.1\t0.02\t5\t6.5\t0\t1\t0\t1\n"
        )

    async def _noop_async(*_args, **_kwargs):
        return None

    original_pl_read_csv = cond_mod.pl.read_csv

    def _patched_pl_read_csv(source, **kwargs):
        if isinstance(source, (str, Path)) and str(source).endswith(".zst"):
            source = io.BytesIO(Path(source).read_bytes())
        return original_pl_read_csv(source, **kwargs)

    monkeypatch.setattr(cond_mod, "meta_analyse_regenie_async", _fake_meta)
    monkeypatch.setattr(cond_mod.pl, "read_csv", _patched_pl_read_csv)
    monkeypatch.setattr(ConditioningStep, "_upload_meta_outputs_to_cloud", _noop_async)
    monkeypatch.setattr(ConditioningStep, "persist_state", _noop_persist)
    asyncio.run(step._meta_analyse_singles())

    summary = pl.read_csv(step.step_path / "summary_output", separator="\t", null_values="NA", infer_schema=False)
    assert summary.select("ID", "CHROM", "LOG10P_OUT").rows() == [("chr1:120:A:G", "1", "6.5")]


def _two_cohort_run(tmp_path, phase):
    """Return a two-cohort run at `phase`."""
    from fedcat.core.engine import DiscoveryRun

    disc = DiscoveryRun(
        config=make_config(local_run_path=str(tmp_path)),
        run_path=tmp_path / "run",
        variant_density_blockfile=tmp_path / "blocks.tsv",
        current_phase=phase,
    )
    disc.platforms = {BioBank.AOU: FakePlatform(), BioBank.UKB: FakePlatform()}
    return disc


@pytest.mark.parametrize(
    "listed, manifest, ok",
    [
        (["out.regenie.gz"], "out.regenie.gz\n", True),
        ([], "", True),  # a block with no annotated variants in the stratum lists nothing
        (["out.regenie.gz"], "out.regenie.gz\nlost.regenie.gz\n", False),
        (None, "", False),  # no manifest at all
    ],
)
def test_step2_outputs_are_checked_against_each_blocks_manifest(tmp_path, listed, manifest, ok):
    from fedcat.core.job import Job
    from fedcat.core.structures import Phase

    disc = _two_cohort_run(tmp_path, Phase.STEP2_AGG)
    job = Job(
        biobank=BioBank.AOU, description="agg", recipe={"op": "step2_agg", "units": [[0, "eur", "1", 1, 100]]}
    )
    disc.jobs = [job]
    (name,) = disc._expected_step2_output_lists(is_singles=False)[BioBank.AOU]
    assert name == "testrun_eur_chr1_1_100_annot1_step2_aggregates.output_filelist"
    folder = f"{disc.config.cohorts[BioBank.AOU].output_location}/meta_analysis/aggregates"
    platform = disc.platforms[BioBank.AOU]
    if listed is not None:
        platform.listing[folder] = [name, *listed]
        platform.files[folder] = {name: manifest}
    assert asyncio.run(disc._verify_step2_output_lists(is_singles=False)) is ok


def test_other_cohorts_results_are_staged_into_the_meta_bucket(tmp_path):
    from fedcat.core.structures import Phase

    disc = _two_cohort_run(tmp_path, Phase.META)
    ukb_root = f"{disc.config.cohorts[BioBank.UKB].output_location}/meta_analysis"
    disc.platforms[BioBank.UKB].files = {
        f"{ukb_root}/single_variants": {"a_TEST_PHENO.regenie.gz": "a", "b_OTHER_PHENO.regenie.gz": "b"},
        f"{ukb_root}/aggregates": {"c_TEST_PHENO.regenie.gz": "c"},
    }
    assert asyncio.run(disc._stage_cohort_results_to_meta())
    meta_base = disc._meta_base_folder_gcs()
    staged = sorted((f.name, dst) for f, dst in disc.platforms[BioBank.AOU].uploaded)
    assert staged == [
        ("a_TEST_PHENO.regenie.gz", f"{meta_base}/ukb_results/single_variants"),
        ("c_TEST_PHENO.regenie.gz", f"{meta_base}/ukb_results/aggregates"),
    ]
    assert disc.platforms[BioBank.UKB].uploaded == []


def _phase_loop(tmp_path, monkeypatch, outcomes: dict, **run_kwargs):
    """Run DiscoveryRun.run over three phases, each phase doing what *outcomes* says; returns (run, phases run, resets)."""
    from fedcat.core.engine import DiscoveryRun
    from fedcat.core.structures import Phase

    phases = [Phase.STEP1, Phase.STEP2_SINGLES, Phase.STEP2_AGG]
    disc = DiscoveryRun(
        config=make_config(local_run_path=str(tmp_path)),
        run_path=tmp_path / "run",
        variant_density_blockfile=tmp_path / "blocks.tsv",
        phases_to_run=phases,
        current_phase=run_kwargs.pop("start", Phase.STEP1),
    )
    disc.had_error = run_kwargs.pop("had_error", False)
    ran: list[tuple[Phase, bool]] = []

    async def _run_phase(self, *, reset_jobs, reset_jobs_list, attempt_job_reset):
        ran.append((self.current_phase, attempt_job_reset))
        match outcomes.get(self.current_phase):
            case "error":
                self.had_error = True
            case "raise":
                raise RuntimeError("phase blew up")
            case "stop":
                self._stop_requested = True

    async def _noop(*_a, **_k):
        return None

    monkeypatch.setattr(DiscoveryRun, "_run_phase", _run_phase)
    monkeypatch.setattr(DiscoveryRun, "persist_state", lambda *_a, **_k: None)
    monkeypatch.setattr(DiscoveryRun, "flush_state", _noop)
    asyncio.run(disc.run(**run_kwargs))
    return disc, [p for p, _ in ran], [r for _, r in ran]


def test_phase_loop_runs_every_phase_in_order(tmp_path, monkeypatch):
    from fedcat.core.structures import Phase

    disc, ran, _ = _phase_loop(tmp_path, monkeypatch, {})
    assert ran == [Phase.STEP1, Phase.STEP2_SINGLES, Phase.STEP2_AGG]
    assert (disc.required, disc.had_error) == (False, False)


@pytest.mark.parametrize(
    "outcome, had_error",
    [("error", True), ("raise", True), ("stop", False)],  # stop: --cond-one-step, resumable without an error
)
def test_phase_loop_stops_where_a_phase_did(tmp_path, monkeypatch, outcome, had_error):
    from fedcat.core.structures import Phase

    disc, ran, _ = _phase_loop(tmp_path, monkeypatch, {Phase.STEP2_SINGLES: outcome})
    assert ran == [Phase.STEP1, Phase.STEP2_SINGLES]
    assert (disc.current_phase, disc.had_error, disc.required) == (Phase.STEP2_SINGLES, had_error, True)


def test_phase_loop_one_phase(tmp_path, monkeypatch):
    from fedcat.core.structures import Phase

    disc, ran, _ = _phase_loop(tmp_path, monkeypatch, {}, one_phase=True)
    assert ran == [Phase.STEP1] and disc.current_phase is Phase.STEP1


def test_phase_loop_retries_a_failed_phase_with_job_reset(tmp_path, monkeypatch):
    from fedcat.core.structures import JobState, Phase

    _, ran, resets = _phase_loop(tmp_path, monkeypatch, {}, start=Phase.STEP2_AGG, had_error=True)
    assert (ran, resets) == ([Phase.STEP2_AGG], [True])
    # without a previous error, a job reset needs --force-reset
    _, ran, _ = _phase_loop(tmp_path, monkeypatch, {}, reset_jobs=JobState.EMPTY)
    assert ran == []
    _, ran, resets = _phase_loop(
        tmp_path, monkeypatch, {}, start=Phase.STEP2_AGG, reset_jobs=JobState.EMPTY, force_reset=True
    )
    assert (ran, resets) == ([Phase.STEP2_AGG], [True])
