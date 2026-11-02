"""Federated forward-selection conditioning of one phenotype.

A [`ConditioningStep`][fedcat.core.conditioning.ConditioningStep] runs one step for
one phenotype and conditioning type: each cohort stratum tests the remaining units
conditional on those selected so far, the results are meta-analysed, and the next
lead unit per chromosome is selected for the following step.
"""

from __future__ import annotations

import asyncio
import gzip
import shutil
import weakref
from pathlib import Path

import polars as pl
from attrs import Factory, define, field
from loguru import logger

from fedcat.core.io import (
    ACAT_OUTPUT_DTYPES,
    COLLINEARITY_DTYPES,
    COLLINEARITY_HEADER,
    META_OUTPUT_DTYPES,
    SUMMARY_OUTPUT_DTYPES,
    TOPHITS_DTYPES,
    extract_tarball_async,
    reformat_collinearity_output_async,
    reformat_mancojo_output_async,
)
from fedcat.core.job import Job
from fedcat.core.logging import NCPU, POLARS_THREADS, Permits, periodic_progress
from fedcat.core.persistence import PersistentState
from fedcat.core.structures import (
    ALL_CHROMS,
    BIOBANK_PROFILES,
    BioBank,
    CondType,
    DiscoveryConfig,
    JobState,
    Phenotype,
    StepError,
    StepState,
    rebase_path,
    stratum_descrip,
)
from fedcat.core.tracking import progress_line
from fedcat.platforms.base import Platform
from fedcat.platforms.gcs_cromwell import gcs
from fedcat.stats.acat import compute_acat
from fedcat.stats.meta_analysis import meta_analyse_regenie_async, p_to_log10p


def _normalise_chr_in_variant_ids() -> pl.Expr:
    """Expression stripping a leading `chr` from `chr:pos:ref:alt` IDs, leaving other IDs unchanged."""
    return pl.col("ID").str.replace(r"(?i)^chr([^:]+:)", "${1}")


def _normalise_chrom_labels() -> pl.Expr:
    """Expression stripping a leading `chr` from the `CHROM` column."""
    return pl.col("CHROM").cast(pl.Utf8).str.replace_all(r"(?i)^chr", "")


def generate_final_output(
    *,
    run_path: Path,
    cond_type: CondType,
    passing_ids: dict[str, list[str]],
    pheno_name: str,
    required: bool,
    end_step: int,
    output_dir: Path | None = None,
    stat_col: str = "LOG10P",
) -> Path | None:
    """Write the conditioning results of a run from its per-step `summary_output` tables.

    Used by [`ConditioningStep.write_final_results`][fedcat.core.conditioning.ConditioningStep.write_final_results]
    and by `fedcat summarise-indep`. Writes a `trace_chr<chrom>_<id>.tsv` per selected
    unit, tracing it across steps, and the summaries `all_sig_conditioned_<type>_<pheno>.tsv`
    (every test of each selected unit, at the last step it was tested) and
    `all_top_sig_conditioned_<type>_<pheno>.tsv` (each selected unit's most significant test).

    Args:
        run_path: Conditioning run directory holding the `<type>_cond_step<N>` directories.
        cond_type: Conditioning type.
        passing_ids: Selected unit IDs per chromosome, in selection order.
        pheno_name: Phenotype name, used in the output filenames.
        required: Whether signals were left, i.e. the results are partial. Partial results
            go to `<type>_cond_partial_results_step<end_step>`, complete ones to
            `<type>_cond_results`.
        end_step: Last step to include.
        output_dir: Results directory to use instead of the default.
        stat_col: Statistic column base name: `"LOG10P"`, or `"P"` for older runs whose
            `summary_output` has `P_IN`/`P_OUT` columns. Kept in the output.

    Returns:
        The results directory, or None if `end_step` is invalid, the results directory
        already exists, or nothing was selected at step 0 (whose `summary_output` is
        then copied to the results directory).

    Raises:
        RuntimeError: If a unit's trace does not match the step it was selected at.
    """
    logger.info("Generating results files")
    if end_step < 0:
        logger.error(f"Invalid end step: {end_step}")
        return None

    stat_in = f"{stat_col}_IN"
    stat_out = f"{stat_col}_OUT"
    # Rename the LOG10P schema keys to the chosen stat column
    read_schema = pl.Schema(
        {
            (stat_in if k == "LOG10P_IN" else stat_out if k == "LOG10P_OUT" else k): v
            for k, v in SUMMARY_OUTPUT_DTYPES.items()
        }
    )

    if output_dir is not None:
        results_dir = output_dir
    elif required:
        results_dir = run_path / f"{cond_type.value}_cond_partial_results_step{end_step}"
    else:
        results_dir = run_path / f"{cond_type.value}_cond_results"

    if results_dir.exists():
        logger.error(f"Results directory {results_dir} already exists, exiting")
        return None
    results_dir.mkdir(parents=True)

    if end_step == 0 and not any(passing_ids.values()):
        step_file = run_path / f"{cond_type.value}_cond_step0" / "summary_output"
        logger.warning(f"No units passed at step0, copying {step_file}")
        shutil.copy(step_file, results_dir / "step0_summary_output")
        return None

    step_files = [run_path / f"{cond_type.value}_cond_step{i}" / "summary_output" for i in range(end_step + 1)]
    passing_id_pairs = [(chrom, id_val) for chrom, ids in passing_ids.items() for id_val in ids]
    filter_df = pl.DataFrame(
        {
            "CHROM": [p[0] for p in passing_id_pairs],
            "GRP_ID": [p[1] for p in passing_id_pairs],
        }
    )

    all_res_dfs = []
    for f in step_files:
        res_df = pl.read_csv(f, separator="\t", schema_overrides=read_schema, null_values="NA")
        match cond_type:
            case CondType.AGGREGATE:
                res_df = res_df.with_columns(
                    pl.col("ID").str.replace_all(r"(\.singleton|\.0\.\d+)$", "").alias("GRP_ID")
                )
            case CondType.SINGLEVAR:
                res_df = res_df.with_columns(pl.col("ID").alias("GRP_ID"))
        res_df = res_df.join(filter_df, on=["CHROM", "GRP_ID"], how="semi")
        all_res_dfs.append(res_df)
    all_res_df = pl.concat(all_res_dfs).sort(["CHROM", "GRP_ID", "STEP"], maintain_order=True)
    del all_res_dfs

    final_df_parts: list[pl.DataFrame] = []
    for chrom, ids in passing_ids.items():
        for i, id in enumerate(ids):
            trace_df = all_res_df.filter((pl.col("CHROM") == chrom) & (pl.col("GRP_ID") == id))
            max_step = trace_df["STEP"].max()
            if max_step != i + 1 and not (i == len(ids) - 1 and max_step == i):
                raise RuntimeError(f"Unit {id} on chr{chrom} passed as #{i + 1} but seen up to step {max_step}")
            trace_df.select(
                [
                    "STEP",
                    "JUMPER",
                    "CHROM",
                    "ID",
                    "TEST",
                    "A1FREQ_IN",
                    "A1FREQ_OUT",
                    "BETA_IN",
                    "BETA_OUT",
                    "SE_IN",
                    "SE_OUT",
                    stat_in,
                    stat_out,
                ]
            ).write_csv(
                results_dir / f"trace_chr{chrom}_{id.replace(':', '-')}.tsv",
                separator="\t",
                null_value="NA",
            )
            jumper_notes = (
                trace_df.filter(pl.col("JUMPER"))
                .with_columns(pl.col("STEP").cast(pl.Utf8))
                .group_by("ID", maintain_order=True)
                .agg(pl.col("STEP").unique(maintain_order=True).alias("_steps"))
                .with_columns((pl.lit("jumped_in_step:") + pl.col("_steps").list.join(";")).alias("NOTE"))
                .select(["ID", "NOTE"])
            )
            final_df = trace_df.filter(pl.col("STEP") == max_step)
            final_df = final_df.join(jumper_notes, on="ID", how="left").with_columns(pl.col("NOTE").fill_null(""))
            final_df = final_df.with_columns(
                [
                    pl.coalesce([pl.col(stat_out), pl.col(stat_in)]).alias(stat_col),
                    pl.coalesce([pl.col("A1FREQ_OUT"), pl.col("A1FREQ_IN")]).alias("A1FREQ"),
                    pl.coalesce([pl.col("BETA_OUT"), pl.col("BETA_IN")]).alias("BETA"),
                    pl.coalesce([pl.col("SE_OUT"), pl.col("SE_IN")]).alias("SE"),
                ]
            )
            final_df = final_df.select(
                ["ID", stat_col, "A1FREQ", "BETA", "SE", "TEST", "CHROM", "START", "END", "STEP", "NOTE"]
            )
            final_df_parts.append(final_df)

    final_df = pl.concat(final_df_parts)
    top_final_df = pl.concat([d.filter(pl.col(stat_col) == d[stat_col].max()) for d in final_df_parts])
    final_df.write_csv(
        results_dir / f"all_sig_conditioned_{cond_type.value}_{pheno_name}.tsv",
        separator="\t",
        null_value="NA",
    )
    top_final_df.write_csv(
        results_dir / f"all_top_sig_conditioned_{cond_type.value}_{pheno_name}.tsv",
        separator="\t",
        null_value="NA",
    )
    if not required:
        logger.success(f"Conditioning complete for {pheno_name} :)")
    return results_dir


# ── ConditioningStep ──────────────────────────────────────────────────────────


@define(kw_only=True)
class ConditioningStep(PersistentState):
    """One conditioning step of one phenotype and conditioning type, as a state machine.

    Created by [`from_discovery_config`][fedcat.core.conditioning.ConditioningStep.from_discovery_config],
    then run with `run_to_completion` and moved to the next step with `increment`. The
    step's state handlers run in `StepState` order: extract the step-0 subset, select
    the next lead units, build and run the cohort jobs, download and check their
    results, and meta-analyse them. State is saved to `state_path` in the step directory.
    """

    config: DiscoveryConfig
    """Run config."""
    base_path: Path
    """`local_run_path`, holding the conditioning run directories."""
    run_prefix: str
    """Run prefix."""
    run_suffix: str
    """`cond_suffix`, naming the conditioning run."""
    cond_type: CondType
    """Conditioning type."""
    input_tarball_path: Path | None = None
    """Step-0 subset tarball; defaults to the phenotype's subset in `base_path`."""
    dx_priority: str = "normal"
    """Default priority of DNAnexus jobs."""
    chromosomes_to_run: list[str] = Factory(lambda: ALL_CHROMS.copy())
    """Chromosomes to condition."""
    pheno: Phenotype
    """Phenotype being conditioned."""
    chroms_cannot_condition: dict[BioBank, dict[str, set[str]]] = Factory(dict)
    """Chromosomes each biobank's strata can no longer test: no units left, or too many
    variants to condition on (`_too_many_variants`)."""
    chroms_ever_conditioned: dict[BioBank, dict[str, set[str]]] = Factory(dict)
    """Chromosomes each biobank's strata have conditioned on at least one unit."""
    step_number: int = 0
    """Step number, from 0."""
    required: bool = True
    """Whether further steps are needed; False once no signals remain."""
    passing_ids: dict[str, list[str]] = Factory(lambda: {chrom: [] for chrom in ALL_CHROMS})
    """Selected unit IDs per chromosome, in selection order."""
    testing_ids: dict[str, list[str]] = Factory(lambda: {chrom: [] for chrom in ALL_CHROMS})
    """Unit IDs per chromosome to test at this step."""
    jobs: list[Job] = Factory(list)
    """Jobs of this step."""
    state: StepState = StepState.READY
    """Current state."""
    prev_state: StepState = StepState.READY
    """State before the last change, restored by `reset_error_state`."""
    n_threads: int | None = None
    """Threads for ACAT; defaults to one fewer than the CPUs."""

    platforms: dict[BioBank, Platform] = field(init=False, factory=dict, repr=False)
    """Platform per biobank, rebuilt from the config on load; not persisted."""

    # ── Creation ──────────────────────────────────────────────────────────────

    @staticmethod
    def _config_derived(config: DiscoveryConfig, pheno_name: str, cond_type: CondType) -> dict:
        """Return the fields a step takes from the config, for creating it or adopting a new config.

        Raises:
            ValueError: If `cond_suffix` is unset.
        """
        if not config.cond_suffix:
            raise ValueError("cond_suffix must be set in the config to run conditioning phases")

        return {
            "config": config,
            "base_path": Path(config.local_run_path),
            "run_prefix": config.run_prefix,
            "run_suffix": config.cond_suffix,
            "chromosomes_to_run": config.chromosomes_to_run,
            "pheno": Phenotype.from_config(config, pheno_name),
        }

    @classmethod
    def from_discovery_config(
        cls,
        config: DiscoveryConfig,
        pheno_name: str,
        cond_type: CondType,
        tarball_path: Path | None = None,
        step_number: int = 0,
    ) -> ConditioningStep:
        """Create a conditioning step for a phenotype.

        Args:
            config: Run config.
            pheno_name: Phenotype to condition.
            cond_type: Conditioning type.
            tarball_path: Step-0 subset tarball; None for the phenotype's subset in
                `local_run_path`.
            step_number: Step number.

        Returns:
            The step.

        Raises:
            ValueError: If `cond_suffix` is unset.
        """
        return cls(
            **cls._config_derived(config, pheno_name, cond_type),
            cond_type=cond_type,
            input_tarball_path=tarball_path,
            step_number=step_number,
        )

    def adopt_config(self, config: DiscoveryConfig) -> None:
        """Replace the persisted config, and the fields derived from it (`--overwrite-config`).

        `input_tarball_path` moves with `base_path` unless it is outside it, and job
        inputs are rebuilt before submission.

        Args:
            config: New config.
        """
        old_base = self.base_path
        for name, value in self._config_derived(config, self.pheno.name, self.cond_type).items():
            setattr(self, name, value)
        self.input_tarball_path = rebase_path(self.input_tarball_path, old_base, self.base_path)
        self.__attrs_post_init__()
        for job in self.jobs:
            job._input_stale = True

    def __attrs_post_init__(self) -> None:
        """Link jobs, build platforms and fill defaults left unset."""
        from fedcat.platforms import build_platforms

        for job in self.jobs:
            job._parent = weakref.ref(self)
        self.platforms = build_platforms(self.config)
        if self.chromosomes_to_run != ALL_CHROMS:
            self.passing_ids = {chrom: self.passing_ids.get(chrom, []) for chrom in self.chromosomes_to_run}
            self.testing_ids = {chrom: self.testing_ids.get(chrom, []) for chrom in self.chromosomes_to_run}
        if self.input_tarball_path is None:
            self.input_tarball_path = (
                self.base_path / f"{self.pheno.name}_{self.cond_type.value}_cond_step0_subset.tar.gz"
            )
        for biobank in self.config.cohorts:
            self.config.cohorts[biobank].output_location = self.config.output_location(biobank)
        if not self.chroms_cannot_condition:
            self.chroms_cannot_condition = {
                bb: {stratum: set() for stratum in cfg.strata} for bb, cfg in self.config.cohorts.items()
            }
        if not self.chroms_ever_conditioned:
            self.chroms_ever_conditioned = {
                bb: {stratum: set() for stratum in cfg.strata} for bb, cfg in self.config.cohorts.items()
            }
        if self.n_threads is None:
            self.n_threads = NCPU - 1

    # ── Path properties ───────────────────────────────────────────────────────

    @property
    def dir_name(self) -> str:
        """Name of the conditioning directory in the cloud outputs; the workflows and applets build the same name."""
        return f"{self.pheno.name}_{self.cond_type.name.lower()}_conditioning_{self.run_suffix}"

    @property
    def run_path(self) -> Path:
        """Local conditioning run directory, prefixed by the run as `local_run_path` may be shared."""
        return self.base_path / f"{self.run_prefix}_{self.dir_name}"

    @property
    def step_path(self) -> Path:
        """Local directory of this step."""
        return self.run_path / f"{self.cond_type.value}_cond_step{self.step_number}"

    @property
    def state_path(self) -> Path:
        """State file of this step."""
        return self.step_path / f"conditioning_state_{self.cond_type.value}_step{self.step_number}.json"

    # ── State serialisation ───────────────────────────────────────────────────

    def _state_file(self) -> Path:
        """Return the state file of this step."""
        return self.state_path

    def _state_queue_key(self) -> str:
        """Return the key this phenotype's state writes are queued under."""
        return f"conditioning_{self.config.run_prefix}_{self.pheno.name}_{self.cond_type.value}"

    def _state_upload_destination(self, path: Path) -> str | None:
        """Return the GCS path to upload a state file to, or None without an AoU cohort."""
        if BioBank.AOU not in self.config.cohorts:
            return None
        return f"{self.config.output_location(BioBank.AOU)}/meta_analysis/{self.dir_name}/{path.name}"

    # ── State management ──────────────────────────────────────────────────────

    def set_state(self, target_state: StepState) -> None:
        """Move to a new state, recording the current one as `prev_state`, and save state.

        Args:
            target_state: State to move to; no change if already in it.
        """
        if target_state != self.state:
            self.prev_state = self.state
            self.state = target_state
            self.persist_state(reason=f"step_state_{self.state.name.lower()}")

    def reset_error_state(
        self,
        target_state: StepState | None = None,
        target_job_state: JobState | None = None,
    ) -> None:
        """Move the step and its failed jobs out of ERROR, or the step into ERROR if it was not in it.

        Args:
            target_state: State to move the step to; None for its state before the error.
            target_job_state: State to move failed jobs to; None for their state before the error.
        """
        if self.state == StepState.ERROR:
            if target_state is not None:
                self.set_state(target_state)
            else:
                self.state = self.prev_state
                self.prev_state = StepState.ERROR
            for job in self.jobs:
                if job.state == JobState.ERROR:
                    job.reset_error_state(target_job_state)
            self.persist_state(reason="error_reset")
        else:
            self.set_state(StepState.ERROR)
            logger.error(f"Attempted to reset ERROR state that was in {self.state.name}")

    # ── Running ───────────────────────────────────────────────────────────────

    async def run(self) -> None:
        """Run the handler of each state in turn until COMPLETED.

        Raises:
            StepError: If the step enters ERROR.
        """
        while self.state != StepState.COMPLETED:
            method_name = f"_{self.state.name.lower()}"
            handler = getattr(self, method_name)
            await handler()

    async def run_to_completion(self) -> bool:
        """Run the step, catching errors.

        Returns:
            True if the step completed or found no further conditioning needed.
        """
        try:
            await self.run()
            return not self.required or self.state == StepState.COMPLETED
        except StepError:
            logger.error("Conditioning step ended in StepError")
            return False
        except (asyncio.CancelledError, KeyboardInterrupt, SystemExit):
            raise
        except BaseException:  # also a fedcat_rs panic, which pyo3 raises as a BaseException
            self.set_state(StepState.ERROR)
            logger.exception("Unexpected error in conditioning")
            return False

    def increment(self) -> None:
        """Select from the completed step's results and move to the next step.

        Per chromosome, results above the retest threshold go forward to the next step,
        whose lead unit is selected in `_sort_input`. A chromosome with no unit below the
        select threshold is finished, and one whose only remaining result is selectable
        has that unit selected outright. Sets `required` to False when nothing goes
        forward. Writes the next step's `tophits_input` before changing this step.

        Raises:
            FileExistsError: If the next step's directory already exists.
        """
        if self.state != StepState.COMPLETED:
            self.set_state(StepState.ERROR)
            logger.error("Cannot increment step that is not in completed state")
            return
        if not self.required:
            self.set_state(StepState.ERROR)
            logger.error("Cannot increment step that is not required")
            return

        logger.info("Processing output")
        step_df = pl.read_csv(
            self.step_path / "summary_output",
            separator="\t",
            schema_overrides=SUMMARY_OUTPUT_DTYPES,
            null_values="NA",
        )

        match self.cond_type:
            case CondType.AGGREGATE:
                step_df = step_df.with_columns(
                    pl.col("ID").str.replace_all(r"(\.singleton|\.0\.\d+)$", "").alias("GRP_ID")
                )
                log10p_retest = p_to_log10p(self.config.cond_p_agg_retest)
                log10p_select = p_to_log10p(self.config.cond_p_agg_select)
            case CondType.SINGLEVAR:
                step_df = step_df.with_columns(pl.col("ID").alias("GRP_ID"))
                log10p_retest = p_to_log10p(self.config.cond_p_singles_retest)
                log10p_select = p_to_log10p(self.config.cond_p_singles_select)
                if self.config.cond_check_coll:
                    coll_df = pl.read_csv(
                        self.step_path / "collinearity" / "all_cohorts_max_collinearity",
                        separator="\t",
                        schema_overrides=COLLINEARITY_DTYPES,
                    )
                    if coll_df.height > 0:
                        step_df = step_df.join(coll_df, on="ID", how="left")
                        step_df = step_df.with_columns(pl.col("max_rsq").fill_null(0))
                        step_df = step_df.filter(pl.col("max_rsq") < self.config.cond_rsq_threshold)

        tophits_df = step_df.filter(pl.col("LOG10P_OUT") > log10p_retest)
        selectable_df = step_df.filter(pl.col("LOG10P_OUT") > log10p_select)

        passed: dict[str, str] = {}
        for aggs in selectable_df.partition_by("CHROM", maintain_order=True):
            chrom = str(aggs["CHROM"][0])
            if tophits_df.filter(pl.col("CHROM") == chrom)["GRP_ID"].n_unique() == 1:
                passed[chrom] = aggs["GRP_ID"][0]
        continuing = [c for c in selectable_df["CHROM"].unique().to_list() if c not in passed]
        tophits_df = tophits_df.filter(pl.col("CHROM").is_in(continuing))

        def _record_passed() -> None:
            """Add the units selected outright to `passing_ids`."""
            for chrom, grp_id in passed.items():
                self.passing_ids[chrom].append(grp_id)
                self.testing_ids[chrom].clear()

        if tophits_df.is_empty():
            _record_passed()
            self.required = False
            return

        next_step_path = self.run_path / f"{self.cond_type.value}_cond_step{self.step_number + 1}"
        if next_step_path.exists():
            logger.error(f"Next step directory {next_step_path} already exists")
            raise FileExistsError(f"The next step directory {next_step_path} already exists")
        next_step_path.mkdir()
        try:
            tophits_df.select(
                ["ID", "LOG10P_OUT", "A1FREQ_OUT", "BETA_OUT", "SE_OUT", "TEST", "CHROM", "START", "END"]
            ).write_csv(next_step_path / "tophits_input", separator="\t", include_header=False, null_value="NA")
        except BaseException:
            shutil.rmtree(next_step_path, ignore_errors=True)
            raise

        logger.info(f"Advancing to step {self.step_number + 1}")
        _record_passed()
        self.step_number += 1
        self.jobs.clear()
        for chrom in self.testing_ids:
            self.testing_ids[chrom].clear()
        self.set_state(StepState.READY)

    # ── Cloud helpers ─────────────────────────────────────────────────────────

    async def _upload_files_to_cloud(self, files: list[Path], destination_subfolder: str) -> None:
        """Copy existing files to the AoU conditioning folder, if there is an AoU cohort; failures are warned."""
        if BioBank.AOU not in self.config.cohorts:
            return
        existing_files = [f for f in files if f.exists() and f.is_file()]
        if not existing_files:
            return
        for biobank in self.config.cohorts:
            if biobank != BioBank.AOU:
                continue
            output_loc = self.config.output_location(biobank)
            destination_folder = f"{output_loc}/meta_analysis/{self.dir_name}/{destination_subfolder}"
            for f in existing_files:
                try:
                    await gcs.upload_file_async(file=f, destination=f"{destination_folder}/{f.name}")
                except Exception as e:  # noqa: BLE001 — storage errors are backend-specific; a failed copy is only warned about
                    logger.warning(f"Failed to upload {f.name} to gs://{destination_folder}: {e}")

    async def _upload_meta_outputs_to_cloud(self) -> None:
        """Copy this step's meta-analysis outputs to the AoU conditioning folder."""
        if self.cond_type == CondType.SINGLEVAR:
            meta_files = [
                self.step_path / "summary_output",
                self.step_path / "META_single_variants" / "META_single_variants.tsv.zst",
            ]
            if self.config.cond_check_coll:
                meta_files.append(self.step_path / "collinearity" / "all_cohorts_max_collinearity")
        else:
            meta_files = [
                self.step_path / "summary_output",
                self.step_path / "acat_output.zst",
                self.step_path / "META_single_variants" / "META_single_variants.tsv.zst",
                self.step_path / "META_single_variants" / "META_single_variants_acat.tsv.zst",
                self.step_path / "META_burden" / "META_burden.tsv.zst",
                self.step_path / "META_mac10" / "META_mac10.tsv.zst",
            ]
        await self._upload_files_to_cloud(
            meta_files, destination_subfolder=f"meta_outputs/step{self.step_number}"
        )

    async def _upload_final_results_to_cloud(self, results_dir: Path) -> bool:
        """Upload the final results to the AoU bucket, returning False if any file failed."""
        if BioBank.AOU not in self.config.cohorts:
            return True
        destination_folder = (
            f"{self.config.output_location(BioBank.AOU)}/meta_analysis/{self.dir_name}/{results_dir.name}"
        )
        ok = True
        for f in (f for f in results_dir.glob("*.tsv") if f.is_file()):
            try:
                await gcs.upload_file_async(file=f, destination=f"{destination_folder}/{f.name}")
            except Exception as e:  # noqa: BLE001 — storage errors are backend-specific; reported through the return value
                logger.error(f"Failed to upload final result {f.name}: {e}")
                ok = False
        return ok

    # ── State handlers, run by `run` for each StepState ──────────────────────

    async def _error(self) -> None:
        """ERROR: stop the step."""
        raise StepError

    async def _ready(self) -> None:
        """READY: at step 0, extract the step-0 subset tarball."""
        if self.step_number == 0:
            logger.info(f"Extracting step0 input: {self.input_tarball_path}")
            dest = self.run_path / f"{self.cond_type.value}_cond_step0"
            if dest.exists():
                logger.warning(f"Step0 directory {dest} already exists, skipping extraction")
            else:
                await extract_tarball_async(self.input_tarball_path, dest)  # type: ignore[arg-type]
        self.set_state(StepState.SORT_INPUT)

    async def _sort_input(self) -> None:
        """SORT_INPUT: select each chromosome's lead unit and set the units to test.

        At step 0 every unit is tested; at later steps the most significant unit is
        selected and the rest are tested conditional on it.
        """
        logger.info("Sorting input")
        if self.step_number == 0:
            with open(self.step_path / "tophits_input") as f:  # noqa: ASYNC230 — reads 10 bytes
                if f.read(10).strip() == "":
                    logger.warning("No tophits in tophits_input, conditioning not required")
                    self.required = False
                    self.set_state(StepState.COMPLETED)
                    return
        tophits = pl.read_csv(
            self.step_path / "tophits_input",
            separator="\t",
            has_header=False,
            schema=TOPHITS_DTYPES,
            null_values="NA",
        )
        if self.step_number == 0:
            tophits = tophits.filter(pl.col("CHROM").is_in(self.chromosomes_to_run))
            tophits.write_csv(
                self.step_path / "tophits_input", separator="\t", include_header=False, null_value="NA"
            )

        match self.cond_type:
            case CondType.AGGREGATE:
                tophits = tophits.with_columns(
                    pl.col("ID").str.replace_all(r"(\.singleton|\.0\.\d+)$", "").alias("GRP_ID")
                )
                log10p_select = p_to_log10p(self.config.cond_p_agg_select)
            case CondType.SINGLEVAR:
                tophits = tophits.with_columns(pl.col("ID").alias("GRP_ID"))
                log10p_select = p_to_log10p(self.config.cond_p_singles_select)
        tophits = tophits.sort("LOG10P_IN", descending=True)

        for hits in tophits.partition_by("CHROM", maintain_order=True):
            chrom = hits["CHROM"][0]
            hits = hits.unique(subset=["GRP_ID"], keep="first", maintain_order=True)
            for hit in hits["GRP_ID"].to_list():
                if hit in self.passing_ids[str(chrom)]:
                    self.set_state(StepState.ERROR)
                    logger.error(f"Unit {hit} in input is also in passing_ids list")
                    return
            if self.step_number == 0:
                self.testing_ids[str(chrom)] = hits["GRP_ID"].to_list()
            else:
                if hits["LOG10P_IN"][0] <= log10p_select:
                    self.set_state(StepState.ERROR)
                    logger.error(f"No unit passing the selection threshold in chr{chrom} input")
                    return
                self.passing_ids[str(chrom)].append(hits["GRP_ID"][0])
                self.testing_ids[str(chrom)] = hits["GRP_ID"][1:].to_list()
        self.set_state(StepState.GENERATE_JOBS)

    async def _generate_jobs(self) -> None:
        """GENERATE_JOBS: upload the IDs to test and build the cohort jobs.

        One WDL job per biobank covering every stratum, or one applet job per stratum
        and chromosome.
        """
        testing_ids_file = self.step_path / f"ids_to_test_step{self.step_number}"
        testing_lines = [
            f"{chrom}\t{grp_id}" for chrom, id_list in self.testing_ids.items() for grp_id in id_list
        ]
        logger.info("Preparing jobs")
        testing_ids_file.write_text("\n".join(testing_lines) + "\n")

        cfg = self.config

        # The returned file reference (DNAnexus file ID or gs:// path) is the job input
        testing_ids_file_descs: dict[BioBank, str] = {}
        for biobank in cfg.cohorts:
            folder = f"{cfg.output_location(biobank)}/meta_analysis/{self.dir_name}/testing_ids_lists"
            ref = await self.platforms[biobank].upload(testing_ids_file, destination=folder, overwrite=True)
            if not ref:
                self.set_state(StepState.ERROR)
                logger.error(f"Could not resolve testing-ids file reference for {biobank.value} cohort")
                return
            testing_ids_file_descs[biobank] = ref

        for bb, cohort in cfg.cohorts.items():
            interface = BIOBANK_PROFILES[bb].workflow_interface

            match interface:
                case "dx_applet_v1":
                    for stratum in cohort.strata:
                        for chrom, test_list in self.testing_ids.items():
                            if test_list and chrom not in self.chroms_cannot_condition[bb][stratum]:
                                job = Job(
                                    biobank=bb,
                                    description=(
                                        f"{stratum_descrip(bb, stratum)}_{self.cond_type.value}"
                                        f"_step{self.step_number}_chr{chrom}"
                                    ),
                                    recipe={
                                        "op": "cond",
                                        "stratum": stratum,
                                        "testing_ids_file_desc": testing_ids_file_descs[bb],
                                        "input_list": [(chrom, ",".join(self.passing_ids[chrom]))],
                                    },
                                    priority=self.dx_priority,
                                )
                                job._parent = weakref.ref(self)
                                self.jobs.append(job)

                case "wdl_scatter_v1":
                    strata = cohort.strata
                    input_list = [
                        (chrom, ",".join(self.passing_ids[chrom]))
                        for chrom, test_list in self.testing_ids.items()
                        if test_list and not all(chrom in self.chroms_cannot_condition[bb][s] for s in strata)
                    ]
                    if input_list:
                        job = Job(
                            biobank=bb,
                            description=(
                                f"{bb.value}_{self.cond_type.value}_step{self.step_number}"
                                f"_chr{','.join(t[0] for t in input_list)}"
                            ),
                            recipe={
                                "op": "cond",
                                "strata": list(strata),
                                "testing_ids_file_desc": testing_ids_file_descs.get(bb, ""),
                                "input_list": input_list,
                                # Each stratum leaves out its own cannot-condition chromosomes
                                "n_tasks": sum(
                                    chrom not in self.chroms_cannot_condition[bb][s]
                                    for s in strata
                                    for chrom, _ in input_list
                                ),
                            },
                        )
                        job._parent = weakref.ref(self)
                        self.jobs.append(job)

                case _:
                    self.set_state(StepState.ERROR)
                    logger.error(f"unsupported workflow_interface {interface!r} for biobank {bb.value!r}")
                    return

        if not self.jobs:
            self.set_state(StepState.ERROR)
            logger.error("No jobs created by _generate_jobs()")
            return
        if len(self.jobs) > 100:
            self.set_state(StepState.ERROR)
            logger.error("Over 100 conditioning jobs created — something wrong")
            return
        self.set_state(StepState.RUN_JOBS)

    async def _job_progress_msg(self) -> str:
        """Return the progress line of this step's jobs, for the periodic progress log."""
        return await progress_line(self.jobs, self.platforms)

    async def _run_jobs(self) -> None:
        """RUN_JOBS: run the jobs, moving on once all have completed."""
        logger.info("Running jobs")
        tasks = [asyncio.create_task(job.run(delay_seconds=float(i))) for i, job in enumerate(self.jobs)]
        async with periodic_progress(self._job_progress_msg):
            await asyncio.gather(*tasks)
        if all(job.state == JobState.COMPLETED for job in self.jobs):
            self.set_state(StepState.DOWNLOAD_RESULTS)
            logger.success("All jobs completed")
        else:
            self.set_state(StepState.ERROR)
            logger.error("Not all jobs completed")

    async def _download_results(self) -> None:
        """DOWNLOAD_RESULTS: download every cohort stratum's results."""
        logger.info("Downloading results")
        tasks = []
        for biobank, cohort_cfg in self.config.cohorts.items():
            for stratum in cohort_cfg.strata:
                tasks.append(self._download_cohort(biobank, stratum))
        await asyncio.gather(*tasks)
        self.set_state(StepState.CHECK_RESULTS)

    async def _download_cohort(self, biobank: BioBank, stratum: str) -> None:
        """Download one stratum's results and marker files for this step."""
        label = stratum_descrip(biobank, stratum)
        logger.debug(f"Beginning download for cohort: {label}")
        output_loc = self.config.output_location(biobank)
        match self.cond_type:
            case CondType.AGGREGATE:
                lovely_regex = (
                    rf"^(.*_{stratum}_chr.*\.regenie.gz|"
                    rf".*_{stratum}_chr.*\.cma\.cojo|"
                    rf".*_{stratum}_chr.*_no_aggregates_remaining.*|"
                    rf".*_{stratum}_chr.*_no_conditioning_variants.*|"
                    rf".*_{stratum}_chr.*_many_variants.*)$"
                )
                folder = f"{output_loc}/meta_analysis/{self.dir_name}/agg_cond_step{self.step_number}/"
            case CondType.SINGLEVAR:
                lovely_regex = (
                    rf"^(.*_{stratum}_chr.*\.regenie.gz|"
                    rf".*_{stratum}_chr.*\.cma\.cojo|"
                    rf".*_{stratum}_chr.*_no_variants_remaining.*|"
                    rf".*_{stratum}_chr.*_no_conditioning_variants.*|"
                    rf".*_{stratum}_chr.*_collinearity.*|"
                    rf".*_{stratum}_chr.*_many_variants.*)$"
                )
                folder = f"{output_loc}/meta_analysis/{self.dir_name}/singles_cond_step{self.step_number}/"
        await self.platforms[biobank].download(
            folder=folder,
            pattern=lovely_regex,
            dest=self.step_path / f"regenie_output_{biobank.value}",
            max_workers=4,
        )
        logger.debug(f"Finished download for cohort: {label}")

    async def _check_results(self) -> None:
        """CHECK_RESULTS: check each stratum returned results for every chromosome it tested."""
        logger.info("Checking results")
        for biobank, cohort_cfg in self.config.cohorts.items():
            for stratum in cohort_cfg.strata:
                chrom_list = [
                    chrom
                    for chrom, ids in self.testing_ids.items()
                    if ids and chrom not in self.chroms_cannot_condition[biobank][stratum]
                ]
                for chrom in chrom_list:
                    match self.cond_type:
                        case CondType.AGGREGATE:
                            self._check_chrom_results_aggs(biobank, stratum, chrom)
                        case CondType.SINGLEVAR:
                            self._check_chrom_results_singles(biobank, stratum, chrom)
        if self.state == StepState.ERROR:
            logger.error("One or more results files had issues")
            return
        self.set_state(StepState.META_ANALYSE)

    def _apply_conditioning_presence(
        self,
        biobank: BioBank,
        stratum: str,
        chrom: str,
        label: str,
        *,
        found_regenie: bool,
        marker_file: Path | None,
    ) -> None:
        """Record whether a stratum conditioned on a chromosome, failing if it stopped doing so.

        Selection is made on the meta-analysis, so a stratum may have no carriers of the
        selected units (e.g. a rare variant driven by a larger stratum). Its job then
        writes a `_no_conditioning_variants` marker and returns unconditioned results,
        which are its conditional results. A stratum that conditioned at an earlier step
        and now reports the marker is an error.
        """
        ever_cond = self.chroms_ever_conditioned[biobank][stratum]
        if marker_file is not None:
            if chrom in ever_cond:
                self.set_state(StepState.ERROR)
                logger.error(
                    f"Chr{chrom} in {label} conditioned at an earlier step but reports no "
                    "conditioning variants at this step"
                )
                return
            try:
                fields = marker_file.read_text().split()
            except OSError:
                fields = []
            cause = fields[0] if fields else "UNKNOWN"
            logger.info(f"Chr{chrom} in {label} ran unconditioned ({cause})")
        elif found_regenie and self.step_number > 0:
            ever_cond.add(chrom)

    def _check_chrom_results_aggs(self, biobank: BioBank, stratum: str, chrom: str) -> None:
        """Check one stratum's aggregate results for a chromosome, by the files and markers returned.

        The job writes `_too_many_variants`, and stops conditioning the chromosome, with
        20,000 or more conditioning variants, fewer than 2 samples per variant, or for a
        binary trait fewer than 10 cases per variant. It writes `_many_variants`, a
        warning, with fewer than 10 samples or 20 cases per variant.
        """
        label = stratum_descrip(biobank, stratum)
        cannot_cond = self.chroms_cannot_condition[biobank][stratum]
        biobank_results_folder = self.step_path / f"regenie_output_{biobank.value}"
        results_files = list(biobank_results_folder.glob(f"*_{stratum}_chr{chrom}_*"))
        if not results_files:
            self.set_state(StepState.ERROR)
            logger.error(f"No results file for {label}_chr{chrom}")
            return
        # An empty result (marker or header-only file) only ends conditioning for the
        # stratum when nothing else came back
        found_empty = False
        found_regenie = False
        no_cond_vars_file: Path | None = None
        for file in results_files:
            if file.name.endswith("_too_many_variants"):
                cannot_cond.add(chrom)
                logger.warning(f"Stopped conditioning on chr{chrom} in {label}")
            elif file.name.endswith("_many_variants"):
                logger.debug(f"Chr{chrom} in {label} has under 10 samples (or 20 cases) per conditioning variant")
            elif file.name.endswith("_no_aggregates_remaining"):
                found_empty = True
            elif file.name.endswith("_no_conditioning_variants"):
                no_cond_vars_file = file
            elif file.name.endswith(".cma.cojo"):
                found_regenie = True
            elif file.name.endswith(".regenie.gz"):
                with gzip.open(file, "rt", encoding="utf-8") as f:
                    line_count = sum(1 for _, _ in zip(f, range(3)))
                if line_count == 2:
                    found_empty = True
                else:
                    found_regenie = True
            else:
                self.set_state(StepState.ERROR)
                logger.error(f"{file} of unknown file type")
        if found_empty:
            if found_regenie:
                logger.debug(f"Chr{chrom} in {label} has no testable aggregates; keeping its other results")
            else:
                cannot_cond.add(chrom)
                logger.debug(f"Chr{chrom} in {label} has no aggregates remaining after conditioning")
        if not found_regenie and chrom not in cannot_cond:
            self.set_state(StepState.ERROR)
            logger.error(f"Chr{chrom} in {label} unexpectedly has no .regenie results")
        self._apply_conditioning_presence(
            biobank, stratum, chrom, label, found_regenie=found_regenie, marker_file=no_cond_vars_file
        )

    def _check_chrom_results_singles(self, biobank: BioBank, stratum: str, chrom: str) -> None:
        """Check one stratum's single-variant results for a chromosome, by the files and markers returned.

        The `_too_many_variants` and `_many_variants` markers are as for aggregates,
        without the case-count limits.
        """
        label = stratum_descrip(biobank, stratum)
        cannot_cond = self.chroms_cannot_condition[biobank][stratum]
        biobank_results_folder = self.step_path / f"regenie_output_{biobank.value}"
        results_files = list(biobank_results_folder.glob(f"*_{stratum}_chr{chrom}_*"))
        if not results_files:
            self.set_state(StepState.ERROR)
            logger.error(f"No results file for {label}_chr{chrom}")
            return
        found_coll = False
        found_regenie = False
        no_cond_vars_file: Path | None = None
        for file in results_files:
            if file.name.endswith("_too_many_variants"):
                cannot_cond.add(chrom)
                logger.warning(f"Stopped conditioning on chr{chrom} in {label}")
            elif file.name.endswith("_many_variants"):
                logger.debug(f"Chr{chrom} in {label} has under 10 samples (or 20 cases) per conditioning variant")
            elif file.name.endswith("_no_variants_remaining"):
                cannot_cond.add(chrom)
                logger.debug(f"Chr{chrom} in {label} has no variants remaining after conditioning")
            elif file.name.endswith("_no_conditioning_variants"):
                no_cond_vars_file = file
            elif file.name.endswith("_collinearity"):
                found_coll = True
            elif file.name.endswith(".cma.cojo"):
                found_regenie = True
            elif file.name.endswith(".regenie.gz"):
                with gzip.open(file, "rt", encoding="utf-8") as f:
                    line_count = sum(1 for _, _ in zip(f, range(2)))
                if line_count == 1:
                    cannot_cond.add(chrom)
                    logger.debug(
                        f"Chr{chrom} in {label} has no variants remaining after conditioning (regenie file empty)"
                    )
                else:
                    found_regenie = True
            else:
                self.set_state(StepState.ERROR)
                logger.error(f"{file} of unknown file type")
        if not found_regenie and chrom not in cannot_cond:
            self.set_state(StepState.ERROR)
            logger.error(f"Chr{chrom} in {label} unexpectedly has no .regenie results")
        if not found_coll and self.config.cond_check_coll and chrom not in cannot_cond:
            self.set_state(StepState.ERROR)
            logger.error(f"Chr{chrom} in {label} unexpectedly has no collinearity file")
        self._apply_conditioning_presence(
            biobank, stratum, chrom, label, found_regenie=found_regenie, marker_file=no_cond_vars_file
        )

    async def _meta_analyse(self) -> None:
        """META_ANALYSE: meta-analyse the step's results and write its `summary_output`."""
        match self.cond_type:
            case CondType.AGGREGATE:
                await self._meta_analyse_aggs()
            case CondType.SINGLEVAR:
                await self._meta_analyse_singles()

    async def _meta_analyse_aggs(self) -> None:
        """Meta-analyse aggregate results, then compute ACAT and compare against the step input.

        Writes `summary_output` with each tested unit's statistics before (`_IN`) and
        after (`_OUT`) this step's conditioning, and marks `JUMPER` units whose -log10 p
        rose by more than 2.
        """
        logger.info("Meta-analysing results")
        for subdir in (
            "META_single_variants",
            "META_burden",
            "META_mac10",
        ):
            subpath = self.step_path / subdir
            if subpath.exists():
                shutil.rmtree(subpath)
            subpath.mkdir()

        burden_cohorts: dict[str, list[Path]] = {}
        singles_cohorts: dict[str, list[Path]] = {}
        mac10_cohorts: dict[str, list[Path]] = {}
        for biobank, cohort_cfg in self.config.cohorts.items():
            regenie_dir = self.step_path / f"regenie_output_{biobank.value}"
            for stratum in cohort_cfg.strata:
                label = stratum_descrip(biobank, stratum)
                burden_glob = list(regenie_dir.glob(f"*_{stratum}_chr*_*_burden_*.regenie.gz"))
                singles_glob = list(regenie_dir.glob(f"*_{stratum}_chr*_singlevars_*.regenie.gz"))
                mac10_glob = list(regenie_dir.glob(f"*_{stratum}_chr*_*_mac10_*.regenie.gz"))
                if self.config.cond_engine == "mancojo" and self.step_number > 0:
                    # From step 1 the mancojo engine writes .cma.cojo files, converted here
                    # to REGENIE's layout. Its .regenie.gz files are unconditional, so are
                    # not used even where a .cma.cojo is missing.
                    burden_glob, singles_glob, mac10_glob = [], [], []
                    for pattern, target in (
                        (f"*_{stratum}_chr*_burden_mancojo.cma.cojo", "burden"),
                        (f"*_{stratum}_chr*_singlevars_mancojo.cma.cojo", "singles"),
                        (f"*_{stratum}_chr*_mac10_mancojo.cma.cojo", "mac10"),
                    ):
                        cma = list(regenie_dir.glob(pattern))
                        if not cma:
                            continue
                        converted = self.step_path / f"mancojo_{label}_{target}.tsv"
                        if not await reformat_mancojo_output_async(cma, converted):
                            continue
                        if target == "burden":
                            burden_glob = [converted]
                        elif target == "singles":
                            singles_glob = [converted]
                        else:
                            mac10_glob = [converted]
                if burden_glob:
                    burden_cohorts[label] = burden_glob
                if singles_glob:
                    singles_cohorts[label] = singles_glob
                if mac10_glob:
                    mac10_cohorts[label] = mac10_glob

        if not burden_cohorts and not singles_cohorts and not mac10_cohorts:
            logger.warning("No aggregates left in any cohort, generating results from input")
            self.required = False
            self.set_state(StepState.COMPLETED)
            return

        meta_tasks = []
        singles_meta = self.step_path / "META_single_variants" / "META_single_variants.tsv.zst"
        # ACAT takes per-bin singles built only from the strata where each variant is
        # rare for the bin, matching the step-0 ACAT built by the meta-analysis job
        acat_singles = singles_meta
        if singles_cohorts:
            meta_tasks.append(meta_analyse_regenie_async(singles_cohorts, singles_meta, "single_variants"))
            if self.config.acat_singles_stratum_rare:
                acat_singles = self.step_path / "META_single_variants" / "META_single_variants_acat.tsv.zst"
                meta_tasks.append(
                    meta_analyse_regenie_async(
                        singles_cohorts, acat_singles, "single_variants", aaf_bins=self.config.aggregate_freqs
                    )
                )
        if burden_cohorts:
            meta_tasks.append(
                meta_analyse_regenie_async(
                    burden_cohorts,
                    self.step_path / "META_burden" / "META_burden.tsv.zst",
                    "burden",
                )
            )
        if mac10_cohorts:
            meta_tasks.append(
                meta_analyse_regenie_async(
                    mac10_cohorts,
                    self.step_path / "META_mac10" / "META_mac10.tsv.zst",
                    "mac10",
                )
            )
        await asyncio.gather(*meta_tasks)

        step0_path = self.run_path / f"{self.cond_type.value}_cond_step0"
        if (step0_path / "annotations").is_file() and (step0_path / "masks").is_file():
            annotations = [step0_path / "annotations"]
            masks = [step0_path / "masks"]
            schema_dir = None
        elif (step0_path / "annotations_parquet").is_dir() and all(
            (step0_path / "annotations_parquet" / f).exists()
            for f in (
                "unit_mask_bounds.parquet",
                "dim_variants.parquet",
                "dim_units.parquet",
                "dim_annots.parquet",
                "mask_annotations.parquet",
                "annotations.parquet",
            )
        ):
            annotations = None
            masks = None
            schema_dir = step0_path / "annotations_parquet"
        else:
            self.set_state(StepState.ERROR)
            logger.error(f"Could not find a complete set of annotations files in {step0_path}")
            return

        # Without single-variant results (e.g. no stratum returned any) ACAT runs on
        # the components that exist
        if not acat_singles.exists():
            acat_singles = self.step_path / "META_single_variants" / "no_singles.tsv"
            acat_singles.write_text("ID\tA1FREQ\tLOG10P\n")
        mac10_meta = self.step_path / "META_mac10" / "META_mac10.tsv.zst"
        burden_meta = self.step_path / "META_burden" / "META_burden.tsv.zst"

        acat_threads = max(self.n_threads, 1)  # type: ignore
        async with Permits(max(acat_threads, POLARS_THREADS, NCPU - 1)):
            await asyncio.to_thread(
                compute_acat,
                single_variants_file=acat_singles,
                burden_file=mac10_meta if mac10_meta.exists() else None,
                schema_dir=schema_dir,
                masks_tsvs=masks,
                annotations_tsvs=annotations,
                output_file=self.step_path / "acat_output.zst",
                aggregate_freqs=self.config.aggregate_freqs,
                n_threads=acat_threads,
                full_bounds=self.config.acat_full_bounds,
            )

        acat_df = pl.read_csv(
            self.step_path / "acat_output.zst",
            separator="\t",
            schema_overrides=ACAT_OUTPUT_DTYPES,
            null_values="NA",
        )
        burden_df = (
            pl.read_csv(
                burden_meta,
                separator="\t",
                comment_prefix="##",
                schema_overrides=META_OUTPUT_DTYPES,
                null_values=["NAN", "NA"],
            )
            if burden_meta.exists()
            else pl.DataFrame(schema=META_OUTPUT_DTYPES)
        )
        tophits_df = pl.read_csv(
            self.step_path / "tophits_input",
            separator="\t",
            has_header=False,
            schema=TOPHITS_DTYPES,
            null_values="NA",
        )

        acat_df = acat_df.rename({"CHR": "CHROM"}).with_columns(
            [
                pl.lit("ACAT").alias("TEST"),
                _normalise_chrom_labels().alias("CHROM"),
                pl.col("LOG10P_ACAT").alias("LOG10P_OUT_ACAT"),
            ]
        )
        burden_df = burden_df.with_columns(
            [
                pl.lit("BURDEN").alias("TEST"),
                _normalise_chrom_labels().alias("CHROM"),
                pl.col("LOG10P").alias("LOG10P_OUT_BURDEN"),
            ]
        ).rename({"A1FREQ": "A1FREQ_OUT", "BETA": "BETA_OUT", "SE": "SE_OUT"})
        tophits_df = tophits_df.with_columns(_normalise_chrom_labels().alias("CHROM"))

        step_df = tophits_df.join(
            acat_df.select(["ID", "TEST", "CHROM", "LOG10P_OUT_ACAT"]), on=["ID", "TEST", "CHROM"], how="left"
        ).join(
            burden_df.select(["ID", "TEST", "CHROM", "A1FREQ_OUT", "BETA_OUT", "SE_OUT", "LOG10P_OUT_BURDEN"]),
            on=["ID", "TEST", "CHROM"],
            how="left",
        )
        step_df = step_df.with_columns(
            [
                pl.coalesce([pl.col("LOG10P_OUT_ACAT"), pl.col("LOG10P_OUT_BURDEN")]).alias("LOG10P_OUT"),
            ]
        ).drop(["LOG10P_OUT_ACAT", "LOG10P_OUT_BURDEN"])
        step_df = step_df.with_columns(
            [
                pl.col("ID").str.replace_all(r"(\.singleton|\.0\.\d+)$", "").alias("AGG"),
                pl.lit(self.step_number).alias("STEP"),
            ]
        )

        # A selected aggregate must not be retested: its ACAT result is dropped, any other is an error
        all_passing = [(chrom, agg) for chrom, ids in self.passing_ids.items() for agg in ids]
        if all_passing:
            passing_df = pl.DataFrame({"CHROM": [p[0] for p in all_passing], "AGG": [p[1] for p in all_passing]})
            flagged = step_df.join(passing_df, on=["CHROM", "AGG"], how="semi").filter(
                pl.col("LOG10P_OUT").is_not_null()
            )
            for row in flagged.filter(pl.col("TEST") == "ACAT").iter_rows(named=True):
                logger.warning(f"Passing agg {row['ID']} found in ACAT results, setting LOG10P_OUT to NA")
            for row in flagged.filter(pl.col("TEST") != "ACAT").iter_rows(named=True):
                self.set_state(StepState.ERROR)
                logger.error(f"Passing agg {row['ID']} found in non-ACAT results (it should not be tested)")
            step_df = (
                step_df.join(
                    passing_df.with_columns(pl.lit(True).alias("_is_passing")), on=["CHROM", "AGG"], how="left"
                )
                .with_columns(
                    pl.when(
                        pl.col("_is_passing") & (pl.col("TEST") == "ACAT") & pl.col("LOG10P_OUT").is_not_null()
                    )
                    .then(pl.lit(None, dtype=pl.Float64))
                    .otherwise(pl.col("LOG10P_OUT"))
                    .alias("LOG10P_OUT")
                )
                .drop("_is_passing")
            )

        step_df = step_df.with_columns((pl.col("LOG10P_OUT") > pl.col("LOG10P_IN") + 2).alias("JUMPER"))
        step_df.drop("AGG").write_csv(self.step_path / "summary_output", separator="\t", null_value="NA")

        for aggs in step_df.partition_by("CHROM", maintain_order=True):
            chrom = aggs["CHROM"][0]
            n_in = aggs.drop_nulls(subset=["LOG10P_IN"]).height
            n_out = aggs.drop_nulls(subset=["LOG10P_OUT"]).height
            if self.step_number == 0:
                if n_out != n_in:
                    logger.warning(f"Missing {n_in - n_out}/{n_in} aggs from meta-analysis for chr{chrom}")
                    logger.debug(
                        f"Missing aggs for chr{chrom}:\n"
                        f"{aggs.filter(pl.col('LOG10P_OUT').is_null()).select(['ID', 'TEST', 'LOG10P_IN', 'LOG10P_OUT'])}"
                    )
            elif n_out <= n_in - 1:
                if n_out == 0:
                    for bb_strata in self.chroms_cannot_condition.values():
                        for cannot_cond in bb_strata.values():
                            cannot_cond.add(str(chrom))
                if n_out < n_in - 1:
                    logger.warning(
                        f"Missing {n_in - 1 - n_out}/{n_in - 1} aggs for chr{chrom} "
                        f"(this is expected if the missing aggs fully overlap the passing aggs)"
                    )
                    logger.debug(
                        f"Passing aggs for chr{chrom}:\n"
                        f"{aggs.filter(pl.col('AGG').is_in(self.passing_ids[str(chrom)])).select(['ID', 'TEST', 'LOG10P_IN', 'LOG10P_OUT'])}"
                    )
                    logger.debug(
                        f"Missing aggs for chr{chrom}:\n"
                        f"{aggs.filter(pl.col('LOG10P_OUT').is_null()).select(['ID', 'TEST', 'LOG10P_IN', 'LOG10P_OUT'])}"
                    )
                if aggs["JUMPER"].any():
                    logger.warning(f"{aggs['JUMPER'].sum()} aggs on chr{chrom} have decreased P-value")
                    logger.debug(
                        f"Jumping aggs for chr{chrom}:\n{aggs.select(['ID', 'TEST', 'LOG10P_IN', 'LOG10P_OUT'])}"
                    )
            else:
                self.set_state(StepState.ERROR)
                logger.error(f"GAINED {1 + n_out - n_in} aggs in meta-analysis for chr{chrom}")

        if self.state == StepState.ERROR:
            return
        await self._upload_meta_outputs_to_cloud()
        self.persist_state(reason=f"meta_outputs_uploaded_step{self.step_number}", max_delay=0)
        self.set_state(StepState.COMPLETED)

    async def _meta_analyse_singles(self) -> None:
        """Meta-analyse single-variant results and collinearity, and compare against the step input.

        Writes `summary_output` with each tested variant's statistics before (`_IN`) and
        after (`_OUT`) this step's conditioning, and marks `JUMPER` variants whose
        -log10 p rose by more than 1.
        """
        logger.info("Meta-analysing results")
        subdirs = ["META_single_variants"]
        if self.config.cond_check_coll:
            subdirs.append("collinearity")
        for subdir in subdirs:
            subpath = self.step_path / subdir
            if subpath.exists():
                shutil.rmtree(subpath)
            subpath.mkdir()

        singles_cohorts: dict[str, list[Path]] = {}
        coll_tasks = []
        for biobank, cohort_cfg in self.config.cohorts.items():
            regenie_dir = self.step_path / f"regenie_output_{biobank.value}"
            for stratum in cohort_cfg.strata:
                label = stratum_descrip(biobank, stratum)
                singles_glob = list(regenie_dir.glob(f"*_{stratum}_chr*_singlevar_*.regenie.gz"))
                if not singles_glob and self.config.cond_engine == "mancojo":
                    # From step 1 the mancojo engine writes .cma.cojo files
                    cma_glob = list(regenie_dir.glob(f"*_{stratum}_chr*_singlevar_*.cma.cojo"))
                    if cma_glob:
                        converted = self.step_path / f"mancojo_{label}_singles.tsv"
                        if await reformat_mancojo_output_async(cma_glob, converted):
                            singles_glob = [converted]
                coll_glob = list(regenie_dir.glob(f"*_{stratum}_chr*_variant_collinearity"))
                coll_path = self.step_path / "collinearity" / f"{label}_singles_collinearity"
                if singles_glob:
                    singles_cohorts[label] = singles_glob
                    if coll_glob:
                        coll_tasks.append(reformat_collinearity_output_async(coll_glob, coll_path))
                    elif self.config.cond_check_coll:
                        coll_path.write_text("\t".join(COLLINEARITY_HEADER) + "\n")

        if not singles_cohorts:
            logger.warning("No variants left in any cohort, generating results from input")
            self.required = False
            self.set_state(StepState.COMPLETED)
            return

        meta_task = meta_analyse_regenie_async(
            singles_cohorts,
            self.step_path / "META_single_variants" / "META_single_variants.tsv.zst",
            "single_variants",
        )
        await asyncio.gather(meta_task, *coll_tasks)

        if self.config.cond_check_coll:
            coll_path = self.step_path / "collinearity"
            coll_files = list(coll_path.glob("*_singles_collinearity"))
            coll_dfs = [pl.read_csv(f, separator="\t", schema_overrides=COLLINEARITY_DTYPES) for f in coll_files]
            coll_df = (
                pl.concat(coll_dfs)
                if coll_dfs
                else pl.DataFrame(schema={"ID": pl.Utf8, "cand_rsq": pl.Float64, "max_rsq": pl.Float64})
            )
            if coll_df.height > 0:
                coll_df = coll_df.group_by("ID").agg(
                    [
                        pl.col("cand_rsq").max(),
                        pl.col("max_rsq").max(),
                    ]
                )
            coll_df.write_csv(coll_path / "all_cohorts_max_collinearity", separator="\t")

        singles_df = pl.read_csv(
            self.step_path / "META_single_variants" / "META_single_variants.tsv.zst",
            separator="\t",
            comment_prefix="##",
            schema_overrides=META_OUTPUT_DTYPES,
            null_values=["NAN", "NA"],
        )
        tophits_df = pl.read_csv(
            self.step_path / "tophits_input",
            separator="\t",
            has_header=False,
            schema=TOPHITS_DTYPES,
            null_values="NA",
        )
        tophits_df = tophits_df.with_columns(_normalise_chrom_labels().alias("CHROM"))

        singles_df = (
            singles_df.with_columns(pl.lit("ADD").alias("TEST"))
            .with_columns(pl.col("ID").str.split(":").list.get(0).alias("CHROM"))
            .with_columns(
                [
                    _normalise_chrom_labels().alias("CHROM"),
                    pl.col("LOG10P").alias("LOG10P_OUT"),
                ]
            )
            .rename({"A1FREQ": "A1FREQ_OUT", "BETA": "BETA_OUT", "SE": "SE_OUT"})
        )
        step_df = tophits_df.join(
            singles_df.select(["ID", "TEST", "CHROM", "A1FREQ_OUT", "BETA_OUT", "SE_OUT", "LOG10P_OUT"]),
            on=["ID", "TEST", "CHROM"],
            how="left",
        ).with_columns(pl.lit(self.step_number).alias("STEP"))

        for vars_df in step_df.partition_by("CHROM", maintain_order=True):
            chrom = vars_df["CHROM"][0]
            output_vars = vars_df.drop_nulls(subset=["LOG10P_OUT"])["ID"].to_list()
            overlap = set(self.passing_ids[str(chrom)]).intersection(set(output_vars))
            if overlap:
                self.set_state(StepState.ERROR)
                logger.error(f"Meta-analysis results present for passing variants {overlap} on chr{chrom}")

        step_df = step_df.with_columns((pl.col("LOG10P_OUT") > pl.col("LOG10P_IN") + 1).alias("JUMPER"))
        step_df.write_csv(self.step_path / "summary_output", separator="\t", null_value="NA")

        for vars_df in step_df.partition_by("CHROM", maintain_order=True):
            chrom = vars_df["CHROM"][0]
            n_in = vars_df.drop_nulls(subset=["LOG10P_IN"]).height
            n_out = vars_df.drop_nulls(subset=["LOG10P_OUT"]).height
            if self.step_number == 0:
                if n_out != n_in:
                    logger.warning(f"Missing {n_in - n_out}/{n_in} vars from meta-analysis for chr{chrom}")
                    logger.debug(
                        f"Missing vars for chr{chrom}:\n"
                        f"{vars_df.filter(pl.col('LOG10P_OUT').is_null()).select(['ID', 'TEST', 'LOG10P_IN', 'LOG10P_OUT'])}"
                    )
            elif n_out <= n_in - 1:
                if n_out == 0:
                    for bb_strata in self.chroms_cannot_condition.values():
                        for cannot_cond in bb_strata.values():
                            cannot_cond.add(str(chrom))
                if n_out < n_in - 1:
                    logger.warning(
                        f"Missing {n_in - 1 - n_out}/{n_in - 1} vars for chr{chrom}"
                        f"(this is expected if missing vars in full ld with passing var)"
                    )
                    logger.debug(
                        f"Passing vars for chr{chrom}:\n"
                        f"{vars_df.filter(pl.col('LOG10P_OUT').is_not_null()).select(['ID', 'TEST', 'LOG10P_IN', 'LOG10P_OUT'])}"
                    )
                    logger.debug(
                        f"Missing vars for chr{chrom}:\n"
                        f"{vars_df.filter(pl.col('LOG10P_OUT').is_null()).select(['ID', 'TEST', 'LOG10P_IN', 'LOG10P_OUT'])}"
                    )
                if vars_df["JUMPER"].any():
                    logger.warning(f"{vars_df['JUMPER'].sum()} vars on chr{chrom} have decreased P-value")
            else:
                self.set_state(StepState.ERROR)
                logger.error(f"GAINED {1 + n_out - n_in} vars in meta-analysis for chr{chrom}")

        if self.state == StepState.ERROR:
            return
        await self._upload_meta_outputs_to_cloud()
        self.persist_state(reason=f"meta_outputs_uploaded_step{self.step_number}", max_delay=0)
        self.set_state(StepState.COMPLETED)

    async def write_final_results(self, partial: bool = False) -> bool:
        """Write the conditioning results and upload them; see [`generate_final_output`][fedcat.core.conditioning.generate_final_output].

        An existing results folder is moved to `<name>_old`.

        Args:
            partial: Whether conditioning stopped with signals left; results go to
                `<type>_cond_partial_results_step<N>`.

        Returns:
            False if the results could not be uploaded.
        """
        # A step that ended before its meta-analysis has no summary, so results end at the one before
        end_step = self.step_number if (self.step_path / "summary_output").exists() else self.step_number - 1
        if end_step < 0:
            return True
        results_dir = self.run_path / (
            f"{self.cond_type.value}_cond_partial_results_step{end_step}"
            if partial
            else f"{self.cond_type.value}_cond_results"
        )
        if results_dir.exists():
            old = results_dir.with_name(f"{results_dir.name}_old")
            logger.warning(f"Moving existing results {results_dir} to {old}")
            if old.exists():
                shutil.rmtree(old)
            results_dir.rename(old)
        built = await asyncio.to_thread(
            generate_final_output,
            run_path=self.run_path,
            cond_type=self.cond_type,
            passing_ids=self.passing_ids,
            pheno_name=self.pheno.name,
            required=partial,
            end_step=end_step,
        )
        return built is None or await self._upload_final_results_to_cloud(built)
