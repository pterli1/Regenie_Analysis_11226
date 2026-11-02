"""[`Job`][fedcat.core.job.Job]: one cloud job of any phase, discovery or conditioning."""

from __future__ import annotations

import asyncio
import weakref
from typing import Any, ClassVar

import dxpy
from attrs import define, field
from loguru import logger

from fedcat.core.structures import (
    BioBank,
    CondType,
    DiscoveryConfig,
    JobKind,
    JobResources,
    JobState,
    stratum_descrip,
    wdl_workflow_name,
    workflow_field,
)
from fedcat.platforms.base import JobSpec, Platform, PlatformError

POLL_INTERVAL = 300
"""Seconds between job status polls."""
RESUME_POLL_DELAY = 0.0
"""Seconds before the first poll of a job resumed in the RUNNING state."""


@define
class Job:
    """One cloud job (an AoU workflow or a UKB applet job) of any phase.

    `recipe` holds just enough to rebuild `input_dict` from the current config, so a
    job picks up config changes made with `--overwrite-config`. Its builder is chosen by
    `recipe["op"]` and the biobank's `workflow_interface`, so a new biobank on an
    existing interface reuses its builders.

    Recipes by `op`, for the `wdl_scatter_v1` (WDL) and `dx_applet_v1` (applet)
    interfaces:

    - `"step1"`: WDL `units` = `[[set, stratum], ...]` in one workflow; applet one job
      per set per stratum, with `stratum` and `set`.
    - `"step2_singles"` / `"step2_agg"`: WDL `units` = `[[group, stratum, chr, start,
      stop], ...]`; applet `stratum`, `block` = `[chr, start, stop, ...]` and `group`.
    - `"meta"`: WDL one workflow over all chromosomes; applet one job per chromosome,
      with `chr`.
    - `"step0_subset"`: `cond_type` (`"singles"` or `"agg"`), `output_folder` and
      `chromosomes`; all chromosomes in one workflow, or one per applet job.
    - `"cond"`: `input_list` = `[(chrom, conditioning_ids), ...]` and
      `testing_ids_file_desc`, the file of IDs to test. WDL one job over all strata
      (`strata`); applet one job per stratum (`stratum`) per chromosome.

    `set` and `group` index `config.step1_pheno_sets` and `config.step2_pheno_groups`.
    An optional `pending` lists the units still to run.
    """

    _parent: weakref.ReferenceType[Any] = field(
        init=False, repr=False, default=None, metadata={"cattrs_exclude": True}
    )
    """Weak reference to the owning `DiscoveryRun` or `ConditioningStep`; not persisted."""

    biobank: BioBank
    """Biobank the job runs on."""
    description: str
    """Unique job name, used in logs, job tags and `--reset-jobs-list`."""
    recipe: dict[str, Any] = field(factory=dict)
    """What the job covers, from which `input_dict` is rebuilt."""
    input_dict: dict[str, Any] = field(factory=dict)
    """Workflow or applet inputs, built by `generate_input`."""
    target: str = ""
    """Workflow name (Cromwell) or applet ID (DNAnexus) to submit."""
    instance_type: str | None = None
    """DNAnexus instance type; None for the applet default."""
    priority: str = "normal"
    """DNAnexus job priority."""
    id: str = ""
    """Platform job ID once submitted."""
    state: JobState = JobState.EMPTY
    """Current state."""
    prev_state: JobState = JobState.EMPTY
    """State before the last change, restored by `reset_error_state`."""

    _input_stale: bool = field(  # type: ignore[assignment]
        init=False, repr=False, default=False, metadata={"cattrs_exclude": True}
    )
    """Whether the config changed since `input_dict` was built; not persisted."""

    @property
    def parent(self) -> Any:
        """The owning `DiscoveryRun` or `ConditioningStep`."""
        return self._parent()  # type: ignore[misc]

    @property
    def config(self) -> DiscoveryConfig:
        """The parent's run config."""
        return self.parent.config

    @property
    def platform(self) -> Platform:
        """The platform of this job's biobank."""
        return self.parent.platforms[self.biobank]

    def set_state(self, target_state: JobState) -> None:
        """Move to a new state, recording the current one as `prev_state`.

        Args:
            target_state: State to move to; no change if already in it.
        """
        if target_state != self.state:
            self.prev_state = self.state
            self.state = target_state

    def reset_error_state(self, target_state: JobState | None = None) -> None:
        """Move a job out of ERROR, or into ERROR if it was not in it.

        Args:
            target_state: State to move to; None returns to the state before the error.
        """
        if self.state == JobState.ERROR:
            if target_state is not None:
                self.set_state(target_state)
            else:
                self.state = self.prev_state
                self.prev_state = JobState.ERROR
        else:
            self.set_state(JobState.ERROR)
            logger.error(
                f"Attempted to reset ERROR state of job {self.description!r}"
                f" that was in {self.prev_state.name} state"
            )

    def _add_covar_stratification(self, prefix: str = "") -> None:
        """Add the covariate-stratification inputs to `input_dict` for a covariate-stratum run.

        Nothing is added for an unstratified run, so it also works with workflow versions
        without these inputs. `prefix` is `"<workflow>."` for WDL inputs and empty for
        applet inputs.
        """
        covar, value = self.config.covar_filter
        if not covar:
            return
        self.input_dict[f"{prefix}stratify_covar"] = covar
        self.input_dict[f"{prefix}stratify_value"] = value

    # Job kind of each recipe op; "cond" depends on the parent step's type
    _KINDS: ClassVar[dict[str, JobKind]] = {
        "step1": JobKind.STEP1,
        "step2_singles": JobKind.STEP2_SINGLES,
        "step2_agg": JobKind.STEP2_AGG,
        "meta": JobKind.META,
        "step0_subset": JobKind.STEP0_SUBSET,
    }

    @property
    def kind(self) -> JobKind | None:
        """What this job runs, or None for a recipe with no known op."""
        op = self.recipe.get("op", "")
        if op == "cond":
            return JobKind.COND_SINGLES if self.parent.cond_type is CondType.SINGLEVAR else JobKind.COND_AGG
        return self._KINDS.get(op)

    @property
    def units(self) -> list[tuple]:
        """Every task this job covers: `(set, stratum)` for step 1, `(group, stratum, chrom, start, stop)`
        for step 2, empty for other ops."""
        r = self.recipe
        op = r.get("op")
        if op not in ("step1", "step2_singles", "step2_agg"):
            return []
        if "units" in r:
            return [tuple(u) for u in r["units"]]
        if op == "step1":
            return [(r.get("set", 0), r["stratum"])]
        chrom, start, stop, *_ = r["block"]
        return [(r.get("group", 0), r["stratum"], str(chrom), str(start), str(stop))]

    @property
    def pending_units(self) -> list[tuple]:
        """The units still to run: `recipe["pending"]` if set, otherwise all `units`."""
        pending = self.recipe.get("pending")
        return self.units if pending is None else [tuple(u) for u in pending]

    @property
    def expected_tasks(self) -> int:
        """Number of platform tasks the job was submitted with, or 0 if unknown.

        The denominator of task progress, as tasks that have not started report no status.
        """
        try:
            if units := self.pending_units:
                return len(units)
        except (AttributeError, KeyError, TypeError):
            return 0  # a job loaded from state whose recipe needs a config to expand
        recipe = self.recipe
        if recipe.get("op") == "step0_subset":
            return len(recipe.get("chromosomes", []))
        if recipe.get("op") == "cond":
            if "n_tasks" in recipe:
                return recipe["n_tasks"]
            return len(recipe.get("input_list", [])) * max(len(recipe.get("strata", [])), 1)
        return 0

    def _grouping_error(self, n_groups: int, used: set[int], what: str) -> bool:
        """Set ERROR and return True if `used` names a group index the config's `n_groups` lacks."""
        if used and max(used) >= n_groups:
            logger.error(
                f"Job {self.description!r} uses {what} {max(used) + 1} but the config now defines {n_groups}; "
                "phenotype grouping changed since the jobs were built, so restart this phase"
            )
            self.set_state(JobState.ERROR)
            return True
        return False

    @property
    def resource_overrides(self) -> dict[str, int | str]:
        """`n_cpu`, `mem_gb` and `disk` set in `input_dict`, to price a job at the size it ran."""
        return {
            name: value
            for name, expected in (("n_cpu", int), ("mem_gb", int), ("disk", str))
            if isinstance(value := self.input_dict.get(f"{self.target}.{name}"), expected)
        }

    def _apply_config_workflow(self) -> None:
        """Set `target`, `instance_type`, `priority` and resource inputs from the current config.

        Applied on every input generation, so a job retried after `--overwrite-config`
        runs on the workflow and machine the new config names. A job whose kind the
        config does not cover is left unchanged.
        """
        kind = self.kind
        if kind is None:
            return
        if workflow := workflow_field(self.config, self.biobank, kind):
            self.target = self.platform.canonical_workflow_ref(workflow)

        resources = self.config.cohorts[self.biobank].workflows.resources.get(kind, JobResources())
        self.instance_type = self.platform.instance_type_for(kind, resources) or self.instance_type
        self.priority = self.platform.priority_for(kind, resources) or self.priority
        self.input_dict.update(self.platform.resource_inputs(self.target, resources))

    async def generate_input(self) -> None:
        """Build `input_dict` and the submission target from the current config.

        Sets the job READY, or ERROR if its inputs cannot be built.
        """
        from fedcat.core.structures import BIOBANK_PROFILES

        op = self.recipe.get("op", "")
        interface = BIOBANK_PROFILES[self.biobank].workflow_interface
        builders = {
            ("step1", "wdl_scatter_v1"): lambda: self._gen_step1_wdl(),
            ("step1", "dx_applet_v1"): lambda: self._gen_step1_dx(),
            ("step2_singles", "wdl_scatter_v1"): lambda: self._gen_step2_wdl(is_singles=True),
            ("step2_agg", "wdl_scatter_v1"): lambda: self._gen_step2_wdl(is_singles=False),
            ("step2_singles", "dx_applet_v1"): lambda: self._gen_step2_dx(is_singles=True),
            ("step2_agg", "dx_applet_v1"): lambda: self._gen_step2_dx(is_singles=False),
            ("meta", "wdl_scatter_v1"): lambda: self._gen_meta_wdl(),
            ("meta", "dx_applet_v1"): lambda: self._gen_meta_dx(),
            ("step0_subset", "wdl_scatter_v1"): lambda: self._gen_subset_wdl(),
            ("step0_subset", "dx_applet_v1"): lambda: self._gen_subset_dx(),
            ("cond", "wdl_scatter_v1"): lambda: self._gen_cond_wdl(),
            ("cond", "dx_applet_v1"): lambda: self._gen_cond_dx(),
        }
        builder = builders.get((op, interface))
        if builder is None:
            logger.error(
                f"No input builder for op={op!r} interface={interface!r} in recipe for {self.description!r}"
            )
            self.set_state(JobState.ERROR)
            return
        builder()
        self._apply_config_workflow()
        self._input_stale = False

    # ── Discovery: STEP 1 ─────────────────────────────────────────────

    def _gen_step1_wdl(self) -> None:
        """Build step 1 WDL inputs for the pending (set, stratum) units."""
        cfg = self.config
        cohort = cfg.cohorts[self.biobank]
        wdl = cohort.workflows.regenie_step1
        wname = wdl_workflow_name(wdl)
        sets, prefixes = cfg.step1_pheno_sets, cfg.step1_run_prefixes
        strata_by_set: dict[int, list[str]] = {}
        for s, stratum in self.pending_units:
            strata_by_set.setdefault(s, []).append(stratum)
        if self._grouping_error(len(sets), set(strata_by_set), "phenotype set"):
            return
        set_ids = sorted(strata_by_set)
        self.input_dict = {
            f"{wname}.phenotype_file": f"gs://{cohort.pheno_file}",
            f"{wname}.covariates_file": f"gs://{cohort.cov_file}",
            f"{wname}.trait_Binary_or_Continuous": cfg.trait_type.value,
            f"{wname}.pheno_transformation_RINT_or_MCC_or_RAW": cfg.trait_transform.value,
            f"{wname}.comma_seperated_list_of_categorical_covariates_or_NONE": cohort.cat_vars,
            f"{wname}.exclcovars": cfg.step1_exclude_covars,
            f"{wname}.pheno_sets": [",".join(sets[s]) for s in set_ids],
            f"{wname}.run_prefixes": [prefixes[s] for s in set_ids],
            f"{wname}.strata_by_set": [strata_by_set[s] for s in set_ids],
            f"{wname}.bucket": cfg.aou_bucket,
            f"{wname}.step1_directory": cohort.step1_directory,
            f"{wname}.step1_input_template": cohort.step1_input_template,
            f"{wname}.regenie_version": cfg.regenie_version,
            f"{wname}.dependencies_path": cfg.dependencies_path(self.biobank),
        }
        self._add_covar_stratification(f"{wname}.")
        self.set_state(JobState.READY)

    def _gen_step1_dx(self) -> None:
        """Build step 1 applet inputs for one phenotype set and stratum."""
        cfg = self.config
        cohort = cfg.cohorts[self.biobank]
        stratum: str = self.recipe["stratum"]
        s: int = self.recipe.get("set", 0)
        if self._grouping_error(len(cfg.step1_pheno_sets), {s}, "phenotype set"):
            return
        self.input_dict = {
            "phenotype_tsv": dxpy.dxlink(cohort.pheno_file),
            "covariate_tsv": dxpy.dxlink(cohort.cov_file),
            "type": cfg.trait_type.value,
            "invnorm": cfg.trait_transform.value,
            "catvars": cohort.cat_vars,
            "exclcovars": cfg.step1_exclude_covars,
            "stratum": stratum,
            "phenos_to_run": ",".join(cfg.step1_pheno_sets[s]),
            "step1_directory": cohort.step1_directory,
            "step1_input_template": cohort.step1_input_template,
            "run_prefix": cfg.step1_run_prefixes[s],
            "regenie_version": cfg.regenie_version,
        }
        self._add_covar_stratification()
        self.set_state(JobState.READY)

    # ── Discovery: STEP 2 ─────────────────────────────────────────────

    def _gen_step2_wdl(self, is_singles: bool) -> None:
        """Build step 2 single-variant or aggregate WDL inputs for the pending units."""
        cfg = self.config
        cohort = cfg.cohorts[self.biobank]
        wdl = cohort.workflows.regenie_step2_singlevar if is_singles else cohort.workflows.regenie_step2_aggregate
        wname = wdl_workflow_name(wdl)
        groups = cfg.step2_pheno_groups
        units_by_group: dict[int, list[list[str]]] = {}
        for g, stratum, chrom, start, stop in self.pending_units:
            units_by_group.setdefault(g, []).append([stratum, chrom, start, stop, cfg.exclude_covars_for(chrom)])
        if self._grouping_error(len(groups), set(units_by_group), "phenotype group"):
            return
        group_ids = sorted(units_by_group)
        cond_file = cohort.cond_file
        self.input_dict = {
            f"{wname}.phenotype_file": f"gs://{cohort.pheno_file}",
            f"{wname}.covariates_file": f"gs://{cohort.cov_file}",
            f"{wname}.trait_Binary_or_Continuous": cfg.trait_type.value,
            f"{wname}.pheno_transformation_RINT_or_MCC_or_RAW": cfg.trait_transform.value,
            f"{wname}.comma_seperated_list_of_categorical_covariates_or_NONE": cohort.cat_vars,
            f"{wname}.pheno_groups": [",".join(groups[g]) for g in group_ids],
            f"{wname}.step1_run_prefixes": [",".join(cfg.step1_prefixes_for(groups[g])) for g in group_ids],
            f"{wname}.group_labels": [cfg.step2_group_label(g) for g in group_ids],
            f"{wname}.units": [units_by_group[g] for g in group_ids],
            f"{wname}.bucket": cfg.aou_bucket,
            f"{wname}.step1_directory": cohort.step1_directory,
            f"{wname}.step1_input_template": cohort.step1_input_template,
            f"{wname}.run_prefix": cfg.run_prefix,
            f"{wname}.results_directory": cohort.root_results_directory,
            f"{wname}.pfile_template": cohort.pfile_template,
            f"{wname}.regenie_version": cfg.regenie_version,
            f"{wname}.dependencies_path": cfg.dependencies_path(self.biobank),
        }
        if cond_file:
            self.input_dict[f"{wname}.conditional_snps_path"] = f"gs://{cond_file}"
        if is_singles:
            self.input_dict[f"{wname}.regenie_mac_threshold"] = str(cfg.singles_mac_threshold)
            if cfg.interaction:
                self.input_dict[f"{wname}.interaction"] = cfg.interaction
        else:
            self.input_dict[f"{wname}.annotation_manifest"] = cohort.annotations_string
            self.input_dict[f"{wname}.coding_only"] = cfg.coding_only
            self.input_dict[f"{wname}.aaf_bins"] = ",".join(cfg.aggregate_freqs)
        self._add_covar_stratification(f"{wname}.")
        self.set_state(JobState.READY)

    def _gen_step2_dx(self, is_singles: bool) -> None:
        """Build step 2 single-variant or aggregate applet inputs for one block."""
        cfg = self.config
        cohort = cfg.cohorts[self.biobank]
        stratum: str = self.recipe["stratum"]
        chrom, start, stop, *_ = self.recipe["block"]
        g: int = self.recipe.get("group", 0)
        if self._grouping_error(len(cfg.step2_pheno_groups), {g}, "phenotype group"):
            return
        phenos = cfg.step2_pheno_groups[g]
        cond_file = cohort.cond_file
        self.input_dict = {
            "phenotype_tsv": dxpy.dxlink(cohort.pheno_file),
            "covariate_tsv": dxpy.dxlink(cohort.cov_file),
            "type": cfg.trait_type.value,
            "invnorm": cfg.trait_transform.value,
            "catvars": cohort.cat_vars,
            "chr": str(chrom),
            "chrom_start": str(start),
            "chrom_end": str(stop),
            "exclcovars": cfg.exclude_covars_for(str(chrom)),
            "stratum": stratum,
            "phenos_to_run": ",".join(phenos),
            "step1_directory": cohort.step1_directory,
            "step1_input_template": cohort.step1_input_template,
            "results_directory": cohort.root_results_directory,
            "run_prefix": cfg.run_prefix,
            "step1_run_prefix": ",".join(cfg.step1_prefixes_for(phenos)),
            "pfile_template": cohort.pfile_template,
            "regenie_version": cfg.regenie_version,
        }
        if label := cfg.step2_group_label(g):
            self.input_dict["pheno_group"] = label
        if cond_file:
            self.input_dict["conditional"] = dxpy.dxlink(cond_file)
        if is_singles:
            self.input_dict["mac"] = str(cfg.singles_mac_threshold)
            if cfg.interaction:
                self.input_dict["interaction"] = cfg.interaction
        else:
            self.input_dict["annotation_manifest"] = cohort.annotations_string
            self.input_dict["coding_only"] = cfg.coding_only
            self.input_dict["aaf_bins"] = ",".join(cfg.aggregate_freqs)
        self._add_covar_stratification()
        self.set_state(JobState.READY)

    # ── Discovery: META ───────────────────────────────────────────────

    def _gen_meta_wdl(self) -> None:
        """Build meta-analysis WDL inputs covering every chromosome."""
        cfg = self.config
        meta_cfg = cfg.meta_analysis_config
        if meta_cfg is None:
            logger.error("meta_analysis_config must be set to build a meta-analysis job")
            self.set_state(JobState.ERROR)
            return
        wdl = meta_cfg.meta_workflow
        wname = wdl_workflow_name(wdl)
        meta_base = f"gs://{self.parent._meta_base_folder_gcs()}"
        input_folders: list[str] = []
        input_biobanks: list[str] = []
        stages = self.parent.meta_mode.stages(cfg.has_cond_file)
        # This job's biobank has its results at the meta base; other cohorts' results are
        # staged under <biobank>_results/ by DiscoveryRun._stage_cohort_results_to_meta
        for bb in cfg.cohorts:
            base = meta_base if bb == self.biobank else f"{meta_base}/{bb.value}_results"
            input_folders.append(f"{base}/single_variants")
            input_biobanks.append(bb.value)
            if stages.aggregates:
                input_folders.append(f"{base}/aggregates")
                input_biobanks.append(bb.value)
        stratum_biobank_pairs = [
            stratum_descrip(bb, stratum) for bb, cohort in cfg.cohorts.items() for stratum in cohort.strata
        ]
        output_folder = f"{meta_base}/all_biobank_meta"

        self.input_dict = {
            f"{wname}.input_folders": input_folders,
            f"{wname}.input_biobanks": input_biobanks,
            f"{wname}.stratum_biobank_pairs": stratum_biobank_pairs,
            f"{wname}.chromosomes": cfg.chromosomes_to_run,
            f"{wname}.output_folder": output_folder,
            f"{wname}.run_prefix": cfg.run_prefix,
            f"{wname}.schema_directory": meta_cfg.schema_directory or "",
            f"{wname}.annotation_manifest": (
                ";".join(str(a) for a in meta_cfg.regenie_annotations) if meta_cfg.regenie_annotations else ""
            ),
            f"{wname}.phenotypes": cfg.pheno_names,
            f"{wname}.create_step0_subset": stages.step0_subset,
            f"{wname}.run_uncond_singles": stages.uncond_singles,
            f"{wname}.run_cond_singles": stages.cond_singles,
            f"{wname}.run_burden": stages.aggregates,
            f"{wname}.run_acat": stages.aggregates,
            f"{wname}.aggs_uncond": not cfg.has_cond_file,
            f"{wname}.acat_full_bounds": cfg.acat_full_bounds,
            f"{wname}.aaf_bins": ",".join(cfg.aggregate_freqs),
            f"{wname}.acat_singles_stratum_rare": cfg.acat_singles_stratum_rare,
            f"{wname}.cond_p_singles_retest": cfg.cond_p_singles_retest,
            f"{wname}.cond_p_agg_retest": cfg.cond_p_agg_retest,
            f"{wname}.step0_singles_folder": f"gs://{self.parent._step0_folder(CondType.SINGLEVAR)}",
            f"{wname}.step0_agg_folder": f"gs://{self.parent._step0_folder(CondType.AGGREGATE)}",
            f"{wname}.dependencies_path": cfg.dependencies_path(self.biobank),
        }
        self.set_state(JobState.READY)

    def _gen_meta_dx(self) -> None:
        """Build meta-analysis applet inputs for one chromosome."""
        cfg = self.config
        meta_cfg = cfg.meta_analysis_config
        if meta_cfg is None:
            logger.error("meta_analysis_config must be set to build a meta-analysis job")
            self.set_state(JobState.ERROR)
            return
        chrom = self.recipe["chr"]
        input_folders: list[str] = []
        input_biobanks: list[str] = []
        stages = self.parent.meta_mode.stages(cfg.has_cond_file)

        for bb, cohort in cfg.cohorts.items():
            input_folders.append(f"{cohort.output_location}/meta_analysis/single_variants")
            input_biobanks.append(bb.value)
            if stages.aggregates:
                input_folders.append(f"{cohort.output_location}/meta_analysis/aggregates")
                input_biobanks.append(bb.value)

        stratum_biobank_pairs = [
            stratum_descrip(bb, stratum) for bb, cohort in cfg.cohorts.items() for stratum in cohort.strata
        ]

        output_folder = f"{cfg.cohorts[self.biobank].output_location}/meta_analysis/all_biobank_meta"

        # The applet takes comma-separated strings, not arrays
        self.input_dict = {
            "input_folders": ",".join(input_folders),
            "input_biobanks": ",".join(input_biobanks),
            "stratum_biobank_pairs": ",".join(stratum_biobank_pairs),
            "chr": str(chrom),
            "output_folder": output_folder.lstrip("/"),
            "run_prefix": cfg.run_prefix,
            "schema_directory": meta_cfg.schema_directory or "",
            "annotation_manifest": (
                ";".join(str(a) for a in meta_cfg.regenie_annotations) if meta_cfg.regenie_annotations else ""
            ),
            "phenotypes": ",".join(cfg.pheno_names),
            "create_step0_subset": stages.step0_subset,
            "run_uncond_singles": stages.uncond_singles,
            "run_cond_singles": stages.cond_singles,
            "run_burden": stages.aggregates,
            "run_acat": stages.aggregates,
            "aggs_uncond": not cfg.has_cond_file,
            "acat_full_bounds": cfg.acat_full_bounds,
            "aaf_bins": ",".join(cfg.aggregate_freqs),
            "acat_singles_stratum_rare": cfg.acat_singles_stratum_rare,
            "cond_p_singles_retest": cfg.cond_p_singles_retest,
            "cond_p_agg_retest": cfg.cond_p_agg_retest,
            "step0_singles_folder": self.parent._step0_folder(CondType.SINGLEVAR).lstrip("/"),
            "step0_agg_folder": self.parent._step0_folder(CondType.AGGREGATE).lstrip("/"),
        }
        self.set_state(JobState.READY)

    # ── Discovery: STEP-0 SUBSET ──────────────────────────────────────

    def _subset_inputs(self) -> dict[str, Any] | None:
        """Return inputs shared by the WDL and applet subset jobs, or None after setting ERROR."""
        cfg = self.config
        meta_cfg = cfg.meta_analysis_config
        if meta_cfg is None:
            logger.error("meta_analysis_config must be set to build a step-0 subset job")
            self.set_state(JobState.ERROR)
            return None
        is_singles = self.recipe["cond_type"] == CondType.SINGLEVAR.value
        # Annotations only subset the aggregate tarball, so a singles job skips downloading them.
        annotations = meta_cfg.regenie_annotations if not is_singles else None
        return {
            "run_prefix": cfg.run_prefix,
            "schema_directory": (meta_cfg.schema_directory or "") if not is_singles else "",
            "annotation_manifest": ";".join(str(a) for a in annotations) if annotations else "",
            "run_uncond_singles": is_singles,
            "run_burden": not is_singles,
            "run_acat": not is_singles,
            "acat_singles_stratum_rare": cfg.acat_singles_stratum_rare,
            "cond_p_singles_retest": cfg.cond_p_singles_retest,
            "cond_p_agg_retest": cfg.cond_p_agg_retest,
        }

    def _gen_subset_wdl(self) -> None:
        """Build step-0 subset WDL inputs for all the recipe's chromosomes."""
        if (shared := self._subset_inputs()) is None:
            return
        wname = wdl_workflow_name(self.config.meta_analysis_config.subset_workflow)
        inputs = shared | {
            "chromosomes": self.recipe["chromosomes"],
            "meta_folder": f"gs://{self.parent._meta_output_folder()}",
            "output_folder": f"gs://{self.recipe['output_folder']}",
            "phenotypes": self.config.pheno_names,
            "dependencies_path": self.config.dependencies_path(self.biobank),
        }
        self.input_dict = {f"{wname}.{k}": v for k, v in inputs.items()}
        self.set_state(JobState.READY)

    def _gen_subset_dx(self) -> None:
        """Build step-0 subset applet inputs for the recipe's single chromosome."""
        if (shared := self._subset_inputs()) is None:
            return
        (chrom,) = self.recipe["chromosomes"]
        self.input_dict = shared | {
            "chr": str(chrom),
            "meta_folder": self.parent._meta_output_folder().lstrip("/"),
            "output_folder": self.recipe["output_folder"].lstrip("/"),
            "phenotypes": ",".join(self.config.pheno_names),
        }
        self.set_state(JobState.READY)

    # ── Conditioning: applet ─────────────────────────────────────────

    def _gen_cond_dx(self) -> None:
        """Build conditioning applet inputs for one stratum and chromosome."""
        step = self.parent  # ConditioningStep
        cfg = self.config
        cohort_cfg = cfg.cohorts[self.biobank]
        stratum: str = self.recipe["stratum"]
        input_list: list[tuple[str, str]] = self.recipe["input_list"]
        if len(input_list) != 1:
            logger.error(f"Each UKB conditioning job must cover exactly one chromosome: {self.description!r}")
            self.set_state(JobState.ERROR)
            return
        chrom = input_list[0][0]
        self.input_dict = {
            "phenotype_tsv": dxpy.dxlink(cohort_cfg.pheno_file),
            "covariate_tsv": dxpy.dxlink(cohort_cfg.cov_file),
            "type": step.pheno.trait_type.value,
            "invnorm": step.pheno.transform.value,
            "catvars": cohort_cfg.cat_vars,
            "chr": chrom,
            "exclcovars": cfg.exclude_covars_for(chrom),
            "stratum": stratum,
            "run_prefix": cfg.run_prefix,
            "step1_run_prefix": cfg.step1_prefixes_for([step.pheno.name])[0],
            "run_suffix": step.run_suffix,
            "pheno_to_run": step.pheno.name,
            "step1_directory": cohort_cfg.step1_directory,
            "step1_input_template": cohort_cfg.step1_input_template,
            "results_directory": cohort_cfg.root_results_directory,
            "pfile_template": cohort_cfg.pfile_template,
            "regenie_version": cfg.regenie_version,
            "cond_step": step.step_number,
        }
        if step.cond_type != CondType.SINGLEVAR and cohort_cfg.cond_file:
            self.input_dict["conditional"] = dxpy.dxlink(cohort_cfg.cond_file)
        match step.cond_type:
            case CondType.AGGREGATE:
                self.input_dict["aggregates"] = dxpy.dxlink(self.recipe["testing_ids_file_desc"])
                self.input_dict["cond_aggregates"] = input_list[0][1]
                self.input_dict["annotation_manifest"] = cohort_cfg.annotations_string
                self.input_dict["aaf_bins"] = ",".join(cfg.aggregate_freqs)
                self.input_dict["min_arm_allele_count"] = str(cfg.cond_min_arm_ac)
                self.input_dict["cond_engine"] = cfg.cond_engine
                self.input_dict["cond_rsq_threshold"] = str(cfg.cond_rsq_threshold)
                self.input_dict["cond_mancojo_effect_mode"] = cfg.cond_mancojo_effect_mode
                self.input_dict["acat_singles_stratum_rare"] = cfg.acat_singles_stratum_rare
            case CondType.SINGLEVAR:
                self.input_dict["vars_to_run"] = dxpy.dxlink(self.recipe["testing_ids_file_desc"])
                self.input_dict["condition_vars"] = input_list[0][1]
                self.input_dict["check_coll"] = str(int(cfg.cond_check_coll))
                self.input_dict["cond_engine"] = cfg.cond_engine
                self.input_dict["cond_rsq_threshold"] = str(cfg.cond_rsq_threshold)
                self.input_dict["cond_mancojo_effect_mode"] = cfg.cond_mancojo_effect_mode
        self._add_covar_stratification()
        self.set_state(JobState.READY)

    # ── Conditioning: WDL ─────────────────────────────────────────────

    def _gen_cond_wdl(self) -> None:
        """Build conditioning WDL inputs covering every stratum and chromosome in one workflow.

        Chromosomes a stratum cannot condition on are left out, as are strata with none left.
        """
        step = self.parent  # ConditioningStep
        cfg = self.config
        cohort_cfg = cfg.cohorts[self.biobank]
        stratum_values: list[str] = self.recipe["strata"]
        input_list: list[tuple[str, str]] = self.recipe["input_list"]
        testing_ids_file_desc: str = self.recipe["testing_ids_file_desc"]

        wdl = (
            cohort_cfg.workflows.conditioning_singlevar
            if step.cond_type == CondType.SINGLEVAR
            else cohort_cfg.workflows.conditioning_aggregate
        )
        wn = wdl_workflow_name(wdl)
        # One list of [chrom, exclcovars, conditioning_ids, testing_ids_file] per stratum
        all_conditioning_blocks: list[list[list[str]]] = []
        strata: list[str] = []
        for stratum in stratum_values:
            cannot_cond = step.chroms_cannot_condition[self.biobank][stratum]
            blocks = [
                [
                    t[0],
                    cfg.exclude_covars_for(t[0]),
                    t[1],
                    testing_ids_file_desc,
                ]
                for t in input_list
                if t[0] not in cannot_cond
            ]
            if blocks:
                all_conditioning_blocks.append(blocks)
                strata.append(stratum)

        if not all_conditioning_blocks:
            logger.error(
                f"No conditioning blocks for AoU job {self.description!r}: every chromosome is excluded "
                "in every stratum"
            )
            self.set_state(JobState.ERROR)
            return

        output_loc = cfg.output_location(self.biobank)
        bucket = f"gs://{output_loc.split('/')[0]}"

        self.input_dict = {
            f"{wn}.phenotype_file": f"gs://{cohort_cfg.pheno_file}",
            f"{wn}.covariates_file": f"gs://{cohort_cfg.cov_file}",
            f"{wn}.trait_Binary_or_Continuous": step.pheno.trait_type.value,
            f"{wn}.pheno_transformation_RINT_or_MCC_or_RAW": step.pheno.transform.value,
            f"{wn}.comma_seperated_list_of_categorical_covariates_or_NONE": cohort_cfg.cat_vars,
            f"{wn}.strata_array": strata,
            f"{wn}.conditioning_blocks": all_conditioning_blocks,
            f"{wn}.bucket": bucket,
            f"{wn}.run_prefix": cfg.run_prefix,
            f"{wn}.step1_run_prefix": cfg.step1_prefixes_for([step.pheno.name])[0],
            f"{wn}.run_suffix": step.run_suffix,
            f"{wn}.pheno_to_run": step.pheno.name,
            f"{wn}.step1_directory": cohort_cfg.step1_directory,
            f"{wn}.step1_input_template": cohort_cfg.step1_input_template,
            f"{wn}.results_directory": cohort_cfg.root_results_directory,
            f"{wn}.pfile_template": cohort_cfg.pfile_template,
            f"{wn}.regenie_version": cfg.regenie_version,
            f"{wn}.dependencies_path": cfg.dependencies_path(self.biobank),
            f"{wn}.cond_step": step.step_number,
        }
        if step.cond_type != CondType.SINGLEVAR and cohort_cfg.cond_file:
            self.input_dict[f"{wn}.conditional_snps_path"] = f"gs://{cohort_cfg.cond_file}"
        match step.cond_type:
            case CondType.AGGREGATE:
                self.input_dict[f"{wn}.annotation_manifest"] = cohort_cfg.annotations_string
                self.input_dict[f"{wn}.aaf_bins"] = ",".join(cfg.aggregate_freqs)
                self.input_dict[f"{wn}.min_arm_allele_count"] = cfg.cond_min_arm_ac
                self.input_dict[f"{wn}.cond_engine"] = cfg.cond_engine
                self.input_dict[f"{wn}.cond_rsq_threshold"] = cfg.cond_rsq_threshold
                self.input_dict[f"{wn}.cond_mancojo_effect_mode"] = cfg.cond_mancojo_effect_mode
                self.input_dict[f"{wn}.acat_singles_stratum_rare"] = cfg.acat_singles_stratum_rare
            case CondType.SINGLEVAR:
                self.input_dict[f"{wn}.check_collinearity"] = cfg.cond_check_coll
                self.input_dict[f"{wn}.cond_engine"] = cfg.cond_engine
                self.input_dict[f"{wn}.cond_rsq_threshold"] = cfg.cond_rsq_threshold
                self.input_dict[f"{wn}.cond_mancojo_effect_mode"] = cfg.cond_mancojo_effect_mode
        self._add_covar_stratification(f"{wn}.")
        self.set_state(JobState.READY)

    async def run(self, delay_seconds: float = 0.0) -> None:
        """Build inputs, submit and wait for the job, from whichever state it is in.

        Returns once the job is COMPLETED or in ERROR. A job already in ERROR is left
        for the caller to reset.

        Args:
            delay_seconds: Seconds to wait before starting.
        """
        if delay_seconds > 0:
            await asyncio.sleep(delay_seconds)
        if self.state == JobState.COMPLETED:
            logger.info(f"Job {self.description!r} already completed, skipping")
            return
        if self.state == JobState.ERROR:
            logger.warning(f"Job {self.description!r} in ERROR state, must be reset first to run")
            return
        resumed = self.state == JobState.RUNNING
        # Each stage sets the next state
        if self.state == JobState.EMPTY:
            await self.generate_input()
        if self.state == JobState.READY:
            await self._submit()
        if self.state == JobState.RUNNING:
            await self._wait(initial_delay=RESUME_POLL_DELAY if resumed else None)

    async def _submit(self) -> None:
        """Submit the job, rebuilding stale inputs first, and persist the parent's state."""
        if self._input_stale or not self.input_dict:
            await self.generate_input()
        if self.state == JobState.ERROR:
            return
        if not self.target:
            self.set_state(JobState.ERROR)
            logger.error(f"Job {self.description!r} has no submission target")
            return
        spec = JobSpec(
            description=self.description,
            target=self.target,
            inputs=self.input_dict,
            scratch_dir=self.parent.run_path,
            instance_type=self.instance_type,
            priority=self.priority,
            tags=[self.description],
        )
        try:
            self.id = await self.platform.submit(spec)
            self.set_state(JobState.RUNNING)
            logger.info(f"Launched {self.biobank.value} job: {self.description!r} ({self.id})")
            self.parent.persist_state(reason=f"launched_{self.description}", max_delay=10)
        except PlatformError as e:
            self.set_state(JobState.ERROR)
            logger.error(f"Failed to launch job {self.description!r}: {e}")
            self.parent.persist_state(reason=f"launch_failed_{self.description}", max_delay=10)

    async def _wait(self, *, initial_delay: float | None = None) -> None:
        """Poll the job until it finishes, then set COMPLETED or ERROR and persist state."""
        status, success = await self.platform.wait(
            self.id, poll_interval=POLL_INTERVAL, description=self.description, initial_delay=initial_delay
        )
        if not success:
            self.set_state(JobState.ERROR)
            logger.error(f"Job {self.description!r} ({self.id}) failed with status {status!r}")
            self.parent.persist_state(reason=f"failed_{self.description}", max_delay=10)
            return
        self.set_state(JobState.COMPLETED)
        logger.success(f"Job {self.description!r} ({self.id}) completed with status {status!r}")
        self.parent.persist_state(reason=f"completed_{self.description}", max_delay=10)
