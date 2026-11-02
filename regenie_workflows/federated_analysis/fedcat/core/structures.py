"""Enums, configuration classes and shared helpers for fedcat.

[`DiscoveryConfig`][fedcat.core.structures.DiscoveryConfig] is the schema of the
config JSON. Runtime-only structures such as
[`Phenotype`][fedcat.core.structures.Phenotype] are derived from it and never persisted.
"""

from __future__ import annotations

import copy
import os
import re
from collections.abc import Mapping, Sequence
from enum import Enum, IntEnum
from pathlib import Path
from typing import Any

from attrs import define, field
from cattrs.preconf.json import make_converter
from loguru import logger

ALL_CHROMS: list[str] = [str(i) for i in range(1, 23)] + ["X"]
MANCOJO_EFFECT_MODES = ("GCTA", "imputeNA", "removeNA")

# Covariate column name, optionally with a baseline level: "VAR" or "VAR[BASE]"
INTERACTION_RE = re.compile(r"[A-Za-z0-9_.\-]+(\[[A-Za-z0-9_.\-<>=+ ]+\])?")
# phenotype names reach the task shell and are comma-joined for --phenoColList
PHENO_NAME_RE = re.compile(r"[A-Za-z0-9_.\-]+")
# stratum names are globbed as _{stratum}_chr and split from their biobank on the first '_'
STRATUM_NAME_RE = re.compile(r"[A-Za-z0-9\-]+")
DX_PRIORITIES: frozenset[str] = frozenset({"high", "normal", "low"})

WORKFLOW_SET_NAME: str = "regenie_wgs_workflow"
DEFAULT_DEPENDENCIES_VERSION: str = "1.7.1"
REGENIE_VERSIONS: tuple[str, ...] = ("4.1.3", "4.1")
DEFAULT_REGENIE_VERSION: str = "4.1.3"

# Covariate stratum label of the unsplit run, and its name on the command line.
COVAR_STRATUM_COMBINED: str = ""
COVAR_STRATUM_COMBINED_NAME: str = "combined"


def compute_dependencies_path(bucket: str, version: str) -> str:
    """Return the `gs://` path of a staged workflow dependencies directory.

    Matches the layout `fedcat setup-workflows` stages and the WDL `dependencies_path`
    defaults: `gs://<bucket>/workflows/<set>/v<version>/dependencies`.

    Args:
        bucket: Workspace bucket name, with or without `gs://`.
        version: Dependencies version, e.g. `"1.7.0"`.

    Returns:
        The dependencies directory path.
    """
    normalised = bucket.removeprefix("gs://").strip("/")
    return f"gs://{normalised}/workflows/{WORKFLOW_SET_NAME}/v{version}/dependencies"


# A workflow name's trailing version, as `fedcat setup-workflows` registers it (-v1-7-0),
# or written with dots
_WORKFLOW_VERSION_SUFFIX = re.compile(r"-v\d+(?:[.-]\d+)*$")


def _workflow_basename(workflow: str) -> str:
    """Strip the directory and `.wdl` extension from a workflow name or path."""
    # Not Path.stem: that would take a version's last ".2" for a file extension
    return workflow.rsplit("/", 1)[-1].removesuffix(".wdl")


def wdl_workflow_name(workflow: str) -> str:
    """Return the name of the workflow inside its WDL, which prefixes the WDL's inputs.

    Args:
        workflow: Workflow name or WDL path, with or without a `-v<version>` suffix.

    Returns:
        `workflow` without directory, `.wdl` extension or version suffix, e.g.
        `regenie_step1_scatter_aou-v1-7-0` gives `regenie_step1_scatter_aou`.
    """
    return _WORKFLOW_VERSION_SUFFIX.sub("", _workflow_basename(workflow))


def versioned_workflow_name(workflow: str, version: str) -> str:
    """Return the Workbench workflow resource name for a workflow at a dependencies version.

    `fedcat setup-workflows` registers each WDL as `<stem>-v<version>` with the
    version's dots as dashes (Workbench rejects dots in workflow ids), so one workspace
    can hold several versions side by side.

    Args:
        workflow: Workflow name or WDL path. A name already ending in a version is
            returned as given.
        version: Dependencies version, e.g. `"1.7.0"`.

    Returns:
        The versioned name, e.g. `regenie_step1_scatter_aou-v1-7-0`.
    """
    name = _workflow_basename(workflow)
    return name if _WORKFLOW_VERSION_SUFFIX.search(name) else f"{name}-v{version.replace('.', '-')}"


def rebase_path(path: Path | None, old_root: Path, new_root: Path) -> Path | None:
    """Move a path from under one root directory to the same place under another.

    Used when `local_run_path` changes: the run's own files follow it, while files the
    user pointed at elsewhere on disk (a tarball or blockfile) are kept.

    Args:
        path: Path to move, or None.
        old_root: Root the path may currently be under.
        new_root: Root to move it to.

    Returns:
        The rebased path, or `path` unchanged when it is None or outside `old_root`.
    """
    if path is None or old_root == new_root:
        return path
    try:
        return new_root / path.relative_to(old_root)
    except ValueError:
        return path


# ── Enum constants ────────────────────────────────────────────────────────────


class BioBank(Enum):
    """Supported biobanks, keyed by their config name."""

    AOU = "aou"
    UKB = "ukb"


class EgressPolicy(Enum):
    """Whether a biobank's result files may be downloaded off its platform.

    `FORBIDDEN` biobanks (e.g. All of Us) may only be run from a process inside that
    biobank's own environment, and their results are never staged to a meta-analysis
    site on another platform.
    """

    FORBIDDEN = "forbidden"
    ALLOWED = "allowed"


@define(frozen=True)
class BiobankProfile:
    """Static properties of a biobank: the platform it runs on and its egress policy.

    Two biobanks with the same `platform_kind` and `workflow_interface` share all
    platform and job-building code and differ only in config.
    """

    platform_kind: str
    """[`Platform.kind`][fedcat.platforms.base.Platform.kind]: `"gcs_cromwell"` or `"dnanexus"`."""
    workflow_interface: str
    """Job-building interface: `"wdl_scatter_v1"` or `"dx_applet_v1"`."""
    egress: EgressPolicy
    """Whether results may leave the platform."""


# Supported biobanks. A biobank reusing an existing platform and interface needs only
# a BioBank member and a row here, plus its cohort config.
BIOBANK_PROFILES: dict[BioBank, BiobankProfile] = {
    BioBank.AOU: BiobankProfile("gcs_cromwell", "wdl_scatter_v1", EgressPolicy.FORBIDDEN),
    BioBank.UKB: BiobankProfile("dnanexus", "dx_applet_v1", EgressPolicy.ALLOWED),
}


class TraitType(Enum):
    """Phenotype trait type (`trait_type` in the config)."""

    BINARY = "Binary"
    CONTINUOUS = "Continuous"


class Transform(Enum):
    """Phenotype transform applied before testing (`trait_transform` in the config)."""

    MCC = "MCC"
    RINT = "RINT"
    RAW = "RAW"


class CondType(Enum):
    """Conditioning analysis type: aggregates or single variants."""

    AGGREGATE = "agg"
    SINGLEVAR = "singles"


class RunEnv(Enum):
    """Execution environment, used for `local_env` and `meta_analysis_config.location`."""

    AOU = "aou"
    LOCAL = "local"


class WorkflowBackend(Enum):
    """Cromwell submission backend for AoU cohorts.

    - `wb_cli`: submit with the Verily Workbench `wb workflow job` CLI.
    - `local_cromwell`: submit to a Cromwell REST server on this machine.
    """

    WB_CLI = "wb_cli"
    LOCAL_CROMWELL = "local_cromwell"


class Phase(Enum):
    """Phases of `fedcat discover`, in run order."""

    STEP1 = "step1"
    STEP2_SINGLES = "step2_singles"
    STEP2_AGG = "step2_aggregates"
    META = "meta_analysis"
    STEP0_SUBSET = "step0_subset"
    COND_SINGLES = "conditioning_singles"
    COND_AGG = "conditioning_aggregates"

    @property
    def label(self) -> str:
        """Short name for the context at the start of log lines."""
        return self._name_.lower()


PHASE_STRINGS_ORDERED: list[str] = [p.value for p in Phase]
PHASE_ORDERED: list[Phase] = list(Phase)


@define(frozen=True)
class MetaStages:
    """Which meta-analysis stages run (the meta workflow's `run_*` inputs)."""

    uncond_singles: bool
    """Unconditional single variants."""
    cond_singles: bool
    """Single variants conditioned on the `cond_file`."""
    aggregates: bool
    """Burden, MAC-10 singletons and ACAT."""
    step0_subset: bool
    """Build the step-0 subsets in the meta job, while its outputs are on disk."""


class MetaMode(Enum):
    """Result types to meta-analyse (`--meta-mode`)."""

    ALL = "all"
    SINGLES_ONLY = "singles_only"
    COND_SINGLES_ONLY = "cond_singles_only"
    UNCOND_SINGLES_ONLY = "uncond_singles_only"
    AGGREGATES_ONLY = "aggregates_only"

    @classmethod
    def from_phases(cls, phases: Sequence[Phase]) -> MetaMode:
        """Pick the mode covering the result types the given phases produce or consume.

        Args:
            phases: Phases being run.

        Returns:
            `SINGLES_ONLY` or `AGGREGATES_ONLY` when the phases only involve one type,
            otherwise `ALL`.
        """
        singles = Phase.STEP2_SINGLES in phases or Phase.COND_SINGLES in phases
        aggregates = Phase.STEP2_AGG in phases or Phase.COND_AGG in phases
        if singles and not aggregates:
            return cls.SINGLES_ONLY
        if aggregates and not singles:
            return cls.AGGREGATES_ONLY
        if not (singles or aggregates) and (Phase.META in phases or Phase.STEP0_SUBSET in phases):
            logger.info(
                "--meta-mode auto: no step2 or conditioning phases specified, so meta-analysing "
                "singles and aggregates by default (pass --meta-mode <mode> to control this)"
            )
        return cls.ALL

    @property
    def runs_aggregates(self) -> bool:
        """Whether this mode meta-analyses aggregates."""
        return self in (MetaMode.ALL, MetaMode.AGGREGATES_ONLY)

    def stages(self, have_cond_file: bool) -> MetaStages:
        """Resolve this mode to the meta-analysis stages to run.

        Conditional singles exist only with a `cond_file`. ACAT uses the conditional
        singles when present, so `AGGREGATES_ONLY` needs the unconditional singles only
        without a `cond_file`. The step-0 subset is built from the unconditional singles
        and the burden/ACAT results.

        Args:
            have_cond_file: Whether the cohorts have a `cond_file`.

        Returns:
            The stages to run.
        """
        return MetaStages(
            uncond_singles=(
                self in (MetaMode.ALL, MetaMode.SINGLES_ONLY, MetaMode.UNCOND_SINGLES_ONLY)
                or (self is MetaMode.AGGREGATES_ONLY and not have_cond_file)
            ),
            cond_singles=have_cond_file and self is not MetaMode.UNCOND_SINGLES_ONLY,
            aggregates=self.runs_aggregates,
            step0_subset=self is not MetaMode.COND_SINGLES_ONLY,
        )


class JobState(IntEnum):
    """State of a cloud job; integer values are used by `--reset-jobs`."""

    ERROR = -1
    EMPTY = 0
    READY = 1
    RUNNING = 2
    COMPLETED = 3


class StepState(IntEnum):
    """State of a conditioning step, in the order its stages run; values used by `--cond-reset-step`."""

    ERROR = -1
    READY = 0
    SORT_INPUT = 1
    GENERATE_JOBS = 2
    RUN_JOBS = 3
    DOWNLOAD_RESULTS = 4
    CHECK_RESULTS = 5
    META_ANALYSE = 6
    COMPLETED = 7


# ── Custom exception types ────────────────────────────────────────────────────


class StepError(Exception):
    """Raised when a conditioning step has entered its ERROR state."""

    def __init__(self, message: str = "Step in ERROR state") -> None:
        """Create the error with a default message."""
        super().__init__(message)


# ── Configuration dataclasses ────────────────────────────────────────────────


@define
class AnnotationEntry:
    """One named REGENIE annotation triplet: mask, annotation and set-list files.

    Each path is the full path in the biobank's cloud storage and may contain a
    `{chr}` placeholder and `.zst` extension, e.g.
    `project-1234:/wgs_annotations/v2026.2/masks_chr{chr}.zst`. No field may contain
    `,` or `;`.
    """

    name: str
    """Annotation name, used to label its aggregate results."""
    mask_file: str
    """REGENIE `--mask-def` file."""
    annot_file: str
    """REGENIE `--anno-file`."""
    setlist_file: str
    """REGENIE `--set-list` file."""

    def __attrs_post_init__(self) -> None:
        """Reject fields containing `,` or `;`, which separate annotation entries."""
        if any(
            sep in attribute
            for sep in [",", ";"]
            for attribute in [self.name, self.mask_file, self.annot_file, self.setlist_file]
        ):
            raise ValueError(f"AnnotationEntry entries cannot contain ',' or ';': {self}")

    def __str__(self) -> str:
        """Comma-separated `name,mask,annot,setlist`, as the workflows take it."""
        return f"{self.name},{self.mask_file},{self.annot_file},{self.setlist_file}"


class JobKind(Enum):
    """What a job runs, and the key its resources are configured under.

    Each value is the name of the config field holding that job's workflow: a field of
    [`CohortWorkflows`][fedcat.core.structures.CohortWorkflows], or of
    [`MetaAnalysisConfig`][fedcat.core.structures.MetaAnalysisConfig] for `META` and
    `STEP0_SUBSET`. `resources` is keyed the same way, and
    [`workflow_field`][fedcat.core.structures.workflow_field] relies on this.
    """

    STEP1 = "regenie_step1"
    STEP2_SINGLES = "regenie_step2_singlevar"
    STEP2_AGG = "regenie_step2_aggregate"
    META = "meta_workflow"
    STEP0_SUBSET = "subset_workflow"
    COND_SINGLES = "conditioning_singlevar"
    COND_AGG = "conditioning_aggregate"


@define
class JobResources:
    """Machine resources for one kind of job; unset fields keep the platform's default.

    The AoU fields are passed as workflow inputs, overriding the WDL defaults. The
    DNAnexus fields apply to applet jobs, whose disk comes with the instance type. Each
    platform ignores the fields that do not apply to it.
    """

    n_cpu: int | None = None
    """AoU task CPUs."""
    mem_gb: int | None = None
    """AoU task memory in GB."""
    disk: str | None = None
    """AoU task disk as `"<GB> SSD"` or `"<GB> HDD"`, e.g. `"200 SSD"`."""
    max_retries: int | None = None
    """AoU retries after a task failure; 0 (the WDL default) fails the job on first error."""
    preemptible: int | None = None
    """AoU attempts on a spot VM before falling back to a standard one; 0 never uses spot."""
    instance_type: str | None = None
    """DNAnexus instance type, e.g. `"mem1_ssd1_v2_x16"`."""
    priority: str | None = None
    """DNAnexus job priority: `"high"` runs on on-demand machines, `"normal"` or `"low"`
    may use spot machines, and the job fails if one is reclaimed."""


def workflow_field(config: DiscoveryConfig, biobank: BioBank, kind: JobKind) -> str:
    """Return the workflow or applet a config names for a kind of job.

    Args:
        config: Run config.
        biobank: Cohort to read the workflow from. Ignored for the meta and step-0
            subset jobs, which are configured once in `meta_analysis_config`.
        kind: Job kind, whose value is the field to read.

    Returns:
        The workflow name or applet ID, or `""` if unset.
    """
    if kind in (JobKind.META, JobKind.STEP0_SUBSET):
        return getattr(config.meta_analysis_config, kind.value, "") if config.meta_analysis_config else ""
    return getattr(config.cohorts[biobank].workflows, kind.value)


@define
class CohortWorkflows:
    """A cohort's workflows, and the resources to run them on.

    AoU values are Workbench workflow names (or WDL names with `local_cromwell`); UKB
    values are DNAnexus applet IDs.
    """

    regenie_step1: str
    """REGENIE step 1 workflow."""
    regenie_step2_singlevar: str
    """REGENIE step 2 single-variant workflow."""
    regenie_step2_aggregate: str
    """REGENIE step 2 aggregate workflow."""
    conditioning_singlevar: str
    """Single-variant conditioning workflow."""
    conditioning_aggregate: str
    """Aggregate conditioning workflow."""
    resources: dict[JobKind, JobResources] = field(factory=dict)
    """Resource overrides per job kind; kinds left out use the platform default
    (`DEFAULT_INSTANCE_TYPES` on DNAnexus, the WDL's own values on AoU)."""


@define
class CohortCfg:
    """Configuration of one biobank's cohort, shared by all its strata.

    Paths are on the biobank's own storage: `gs://` paths (prefix optional) for AoU,
    DNAnexus file IDs or `project-xxx:/path` for UKB.
    """

    workflows: CohortWorkflows
    """Workflows to run, and their resources."""
    pheno_file: str
    """Phenotype file holding every stratum's individuals."""
    cov_file: str
    """Covariate file."""
    strata: list[str]
    """Sample strata to run and meta-analyse (e.g. genetic ancestry groups). Each
    stratum's samples are those of its step1 input; names may only contain letters,
    digits and `-`."""
    cat_vars: str
    """Comma-separated categorical covariates, or `NONE`."""
    step1_directory: str
    """Step 1 output directory, relative to the bucket or project root."""
    root_results_directory: str
    """Root directory of all other results, relative to the bucket or project root."""
    pfile_template: str
    """PLINK2 pfile set per chromosome, without extension and with a `{chr}` placeholder."""
    step1_input_template: str
    """Step 1 input file set per stratum, without extension and with a `{stratum}`
    placeholder."""
    regenie_annotations: list[AnnotationEntry] = field(factory=list)
    """Annotation triplets for aggregate testing; only needed for aggregate phases."""
    cond_file: str = ""
    """File of variants to always condition on. Set for all cohorts or none."""
    output_location: str = ""
    """Results directory of this run; set at runtime from
    [`DiscoveryConfig.output_location`][fedcat.core.structures.DiscoveryConfig.output_location],
    overriding any config value."""
    cohort_sizes: dict[str, int] = field(factory=dict)
    """Individuals per stratum, used to size step 2 blocks; missing strata default to
    100,000."""
    workflow_backend: WorkflowBackend = WorkflowBackend.WB_CLI
    """AoU workflow submission backend."""
    dependencies_version: str = DEFAULT_DEPENDENCIES_VERSION
    """Version of the workflow dependencies staged by `fedcat setup-workflows` to run."""
    dx_project_id: str = ""
    """DNAnexus project jobs run and write to; empty uses the top-level `dx_project_id`."""
    applet_project_id: str = ""
    """DNAnexus project holding the applets; empty uses the shared applet project."""
    max_tasks_per_workflow: int = 0
    """Split each step 1/step 2 phase over workflows of at most this many tasks (WDL
    cohorts only); 0 for one workflow."""
    max_concurrent_jobs: int = 0
    """Most jobs running at once; 0 for no limit."""

    def __attrs_post_init__(self) -> None:
        """Check the path templates have their placeholders."""
        if "{chr}" not in self.pfile_template:
            raise ValueError(f"pfile_template must contain '{{chr}}' placeholder: {self.pfile_template!r}")
        if "{stratum}" not in self.step1_input_template:
            raise ValueError(
                f"step1_input_template must contain '{{stratum}}' placeholder: {self.step1_input_template!r}"
            )

    @property
    def annotations_string(self) -> str:
        """Annotation triplets as the workflows take them: `name,mask,annot,set;name,...`."""
        return ";".join(str(a) for a in self.regenie_annotations)


@define
class MetaAnalysisConfig:
    """Meta-analysis configuration, required by the `meta_analysis` and `step0_subset` phases.

    Annotations (one of `schema_directory` or `regenie_annotations`) are only needed
    when the meta-analysis runs its aggregate stages.
    """

    location: RunEnv
    """Where the meta-analysis and step-0 subset jobs run: `aou` (Cromwell and GCS) or
    `local` (DNAnexus). May differ from `local_env`, where fedcat itself runs."""
    meta_workflow: str
    """Meta-analysis workflow name, or applet ID for `local`."""
    subset_workflow: str = ""
    """Step-0 subset workflow name, or applet ID for `local`."""
    regenie_annotations: list[AnnotationEntry] | None = None
    """REGENIE-style annotation triplets for the aggregate stages."""
    schema_directory: str | None = None
    """Directory of annotations converted by `fedcat convert-annots`, with a `{chr}`
    placeholder."""


@define
class CovarStratification:
    """Run the analysis separately within each level of a covariate (`stratify_by_covar`)."""

    covar_name: str
    """Column of each cohort's `cov_file` to stratify on."""
    covar_values_names: dict[str, str]
    """Covariate stratum name for each value to analyse; values are compared as strings
    exactly as they appear in the file."""
    also_run_combined: bool = True
    """Also run the unstratified analysis."""

    def value_for(self, covar_stratum: str) -> str:
        """Return the `covar_name` value selecting a covariate stratum.

        Args:
            covar_stratum: Covariate stratum name, or `""` for the combined run.

        Returns:
            The covariate value, or `""` for the combined run.

        Raises:
            ValueError: If `covar_stratum` is not configured.
        """
        if covar_stratum == COVAR_STRATUM_COMBINED:
            return ""
        for value, name in self.covar_values_names.items():
            if name == covar_stratum:
                return value
        raise ValueError(f"Unknown covariate stratum: {covar_stratum}")


@define
class DiscoveryConfig:
    """Configuration of a federated discovery run: the schema of the config JSON.

    Only the fields a run uses are required: e.g. `meta_analysis_config` only for the
    `meta_analysis` and `step0_subset` phases, `cond_suffix` only for conditioning.
    """

    cohorts: dict[BioBank, CohortCfg]
    """Configuration of each biobank to run."""
    pheno_names: list[str]
    """Phenotypes to run; each phenotype and covariate file needs a column per name."""
    trait_type: TraitType
    """Phenotype trait type."""
    trait_transform: Transform
    """Phenotype transform; binary traits should use `RAW`."""
    chromosomes_to_run: list[str]
    """Chromosomes to test."""
    local_run_path: str
    """Local directory holding run directories, state files, logs and conditioning results."""
    run_prefix: str
    """Name of the run, used in its directory and output names."""
    local_env: RunEnv
    """Environment fedcat runs in. Cohorts that forbid data egress (AoU) need `aou`."""
    meta_analysis_config: MetaAnalysisConfig | None = None
    """Meta-analysis configuration; required by the `meta_analysis` and `step0_subset` phases."""
    stratify_by_covar: CovarStratification | None = None
    """Run one independent analysis per level of a covariate; see
    [`covar_stratum_configs`][fedcat.core.structures.DiscoveryConfig.covar_stratum_configs]."""
    covar_stratum: str = ""
    """Covariate stratum this config runs; set by `covar_stratum_configs`, not in the config file."""
    coding_only: bool = False
    """Restrict aggregate tests to sets whose names start with `coding_`."""
    interaction: str = ""
    """Covariate for GxE interaction testing (REGENIE `--interaction`), as `VAR` or
    `VAR[BASE]` to set a categorical baseline. Step 2 single variants only."""
    step1_run_prefix: str | None = None
    """Prefix of step 1 outputs, to share a step 1 fit between runs; defaults to `run_prefix`."""
    step1_pheno_group_size: int = 0
    """Phenotypes fitted together per step 1 job; 0 fits all together."""
    step2_pheno_group_size: int = 0
    """Phenotypes tested together per step 2 task; 0 tests all together."""
    max_step2_singles_block_length: int | None = None
    """Ceiling on the step 2 single-variant block length in bp; None for the platform default."""
    max_step2_aggregates_block_length: int | None = None
    """Ceiling on the step 2 aggregate block length in bp; None for the platform default."""
    regenie_version: str = DEFAULT_REGENIE_VERSION
    """REGENIE version every REGENIE job runs, one of `REGENIE_VERSIONS`."""
    singles_mac_threshold: int = 11
    """Minimum minor allele count for testing single variants."""
    aggregate_freqs: list[str] = field(factory=lambda: ["0.001"])
    """Upper allele-frequency bounds of aggregate tests, passed verbatim to REGENIE
    `--aaf-bins` and used as the ACAT bin label."""
    acat_full_bounds: bool = True
    """Take ACAT mask START/END from all annotated variants, not only those with
    single-variant results, so burden-only masks get positions."""
    acat_singles_stratum_rare: bool = True
    """Include a stratum's single-variant result in a bin's ACAT only where the variant is
    rare (A1FREQ <= bin) in that stratum, rather than filtering on the meta A1FREQ."""
    cond_suffix: str | None = None
    """Suffix of conditioning run directory names; required by the conditioning phases."""
    pgs_prefix: str | None = None
    """LOCO PGS covariate prefix: column `{pgs_prefix}{chrom}` is excluded on each chromosome."""
    chromosomes_no_pgs: list[str] = field(factory=list)
    """Chromosomes without a PGS covariate to exclude."""
    exclude_covars: list[str] = field(factory=list)
    """Covariates excluded on every chromosome."""
    exclude_covars_dict: dict[str, list[str]] = field(factory=dict)
    """Further covariates excluded per chromosome."""
    dx_project_id: str = "project-GbZ9g7QJVJzZVQj45j6xgxZ8"
    """DNAnexus project UKB jobs run and write to."""
    aou_bucket: str = ""
    """GCS bucket AoU results are written to; chosen interactively and saved if unset."""
    aou_bucket_resource_id: str = ""
    """Workbench resource id of `aou_bucket` (for `wb --bucket-id`); resolved automatically."""
    cond_p_agg_retest: float = 1e-6
    """Aggregate results with p above this are dropped at each conditioning step and
    left out of the step-0 subset."""
    cond_p_agg_select: float = 1e-7
    """Aggregates with p below this can be selected; must be <= `cond_p_agg_retest`."""
    cond_p_singles_retest: float = 1e-6
    """Single-variant results with p above this are dropped at each conditioning step
    and left out of the step-0 subset."""
    cond_p_singles_select: float = 5e-8
    """Single variants with p below this can be selected; must be <= `cond_p_singles_retest`."""
    cond_rsq_threshold: float = 0.9
    """Most r² a single variant may have with the selected set to be selected; also the
    COJO collinearity cutoff for the `mancojo` engine."""
    cond_check_coll: bool = True
    """Apply the `cond_rsq_threshold` collinearity check to single variants."""
    cond_max_steps: int = -1
    """Most conditioning steps; -1 for no limit."""
    cond_min_arm_ac: int = 0
    """Binary traits on the `regenie` engine: minor-allele copies required in each of
    cases and controls for a conditioning variant to be kept as a covariate."""
    cond_engine: str = "regenie"
    """How conditional statistics are derived: `regenie` reruns REGENIE with the selected
    variants as covariates, `mancojo` projects marginal statistics through an LD matrix."""
    cond_mancojo_effect_mode: str = "imputeNA"
    """manc-cojo `--effect-mode` for the `mancojo` engine, one of `MANCOJO_EFFECT_MODES`."""

    def __attrs_post_init__(self) -> None:
        """Check the conditioning settings and `cond_file`s, and default `aou_bucket`."""
        if self.cond_engine not in ("regenie", "mancojo"):
            raise ValueError(f"cond_engine must be 'regenie' or 'mancojo', got {self.cond_engine!r}")
        if self.cond_mancojo_effect_mode not in MANCOJO_EFFECT_MODES:
            raise ValueError(
                f"cond_mancojo_effect_mode must be one of {MANCOJO_EFFECT_MODES}, got {self.cond_mancojo_effect_mode!r}"
            )
        for kind in ("agg", "singles"):
            retest, select = getattr(self, f"cond_p_{kind}_retest"), getattr(self, f"cond_p_{kind}_select")
            if select > retest:
                raise ValueError(
                    f"cond_p_{kind}_select ({select:g}) must not exceed cond_p_{kind}_retest ({retest:g})"
                )
        if not self.aou_bucket:
            bucket = os.environ.get("WORKSPACE_BUCKET")
            if bucket:
                self.aou_bucket = bucket

        cond_files = [c.cond_file for c in self.cohorts.values()]
        if any(cond_files) and not all(cond_files):
            raise ValueError("If cond_file is provided for one cohort, it must be provided for all cohorts")

    @property
    def has_cond_file(self) -> bool:
        """Whether the cohorts have a `cond_file` (set for all cohorts or none)."""
        return all(c.cond_file for c in self.cohorts.values())

    @property
    def discovery_run_path(self) -> Path:
        """Local directory holding this run's state files, logs and downloads."""
        return Path(self.local_run_path) / f"{self.run_prefix}_discovery_run"

    @property
    def effective_step1_run_prefix(self) -> str:
        """Prefix of step 1 outputs: `step1_run_prefix`, or `run_prefix` if unset."""
        return self.step1_run_prefix or self.run_prefix

    def _pheno_groups(self, size: int) -> list[list[str]]:
        """Split `pheno_names` into groups of `size`; one group if `size <= 0`."""
        if size <= 0:
            return [list(self.pheno_names)]
        return [self.pheno_names[i : i + size] for i in range(0, len(self.pheno_names), size)]

    @property
    def step1_pheno_sets(self) -> list[list[str]]:
        """Phenotypes fitted together in each step 1 job, in `pheno_names` order."""
        return self._pheno_groups(self.step1_pheno_group_size)

    @property
    def step1_run_prefixes(self) -> list[str]:
        """Output prefix of each step 1 set; a single set keeps the bare prefix."""
        sets = self.step1_pheno_sets
        if len(sets) == 1:
            return [self.effective_step1_run_prefix]
        return [f"{self.effective_step1_run_prefix}_set{i + 1}" for i in range(len(sets))]

    def step1_prefixes_for(self, phenos: Sequence[str]) -> list[str]:
        """Return the prefixes of the step 1 sets holding any of the given phenotypes.

        Args:
            phenos: Phenotype names.

        Returns:
            Step 1 output prefixes, in set order.
        """
        wanted = set(phenos)
        return [
            prefix
            for prefix, pheno_set in zip(self.step1_run_prefixes, self.step1_pheno_sets)
            if wanted.intersection(pheno_set)
        ]

    @property
    def step2_pheno_groups(self) -> list[list[str]]:
        """Phenotypes tested together in each step 2 task, in `pheno_names` order."""
        return self._pheno_groups(self.step2_pheno_group_size)

    def step2_group_label(self, group: int) -> str:
        """Return the tag naming a step 2 phenotype group's outputs.

        Args:
            group: Index into `step2_pheno_groups`.

        Returns:
            `g<N>` (1-based), or `""` when there is a single group.
        """
        return "" if len(self.step2_pheno_groups) == 1 else f"g{group + 1}"

    @property
    def covar_filter(self) -> tuple[str, str]:
        """`(column, value)` the workflows subset samples to.

        `("", "")` when this config is not a covariate stratum, telling job builders to
        omit the stratification inputs.
        """
        if self.stratify_by_covar is None or self.covar_stratum == COVAR_STRATUM_COMBINED:
            return ("", "")
        return (self.stratify_by_covar.covar_name, self.stratify_by_covar.value_for(self.covar_stratum))

    def covar_stratum_configs(self) -> list[tuple[str, DiscoveryConfig]]:
        """Expand this config into one config per covariate stratum to run.

        Each covariate stratum gets a deep copy with:

        - `run_prefix` and `step1_run_prefix` suffixed with the stratum name, which
          separates its step 1 fit, outputs, run directory and state files;
        - the stratifying covariate added to `exclude_covars`, as it is constant
          within the stratum;
        - the stratifying covariate removed from each cohort's `cat_vars`.

        The combined run, when requested, is an unchanged copy. Call this on the parsed
        config, never on a config it returned, so the result is the same on every resume.

        Returns:
            `(covar_stratum, config)` pairs, the combined run first with stratum `""`.
            `[("", self)]` when `stratify_by_covar` is unset.
        """
        strat = self.stratify_by_covar
        if strat is None:
            return [(COVAR_STRATUM_COMBINED, self)]

        variants: list[tuple[str, DiscoveryConfig]] = []
        if strat.also_run_combined:
            combined = copy.deepcopy(self)
            combined.covar_stratum = COVAR_STRATUM_COMBINED
            variants.append((COVAR_STRATUM_COMBINED, combined))

        for name in strat.covar_values_names.values():
            cfg = copy.deepcopy(self)
            cfg.covar_stratum = name
            cfg.step1_run_prefix = f"{self.effective_step1_run_prefix}_{name}"
            cfg.run_prefix = f"{self.run_prefix}_{name}"
            if strat.covar_name not in cfg.exclude_covars:
                cfg.exclude_covars = [*cfg.exclude_covars, strat.covar_name]
            for cohort in cfg.cohorts.values():
                cohort.cat_vars = _drop_covar(cohort.cat_vars, strat.covar_name)
            variants.append((name, cfg))
        return variants

    def exclude_covars_for(self, chrom: str) -> str:
        """Return REGENIE's `--covarExcludeList` for a chromosome.

        See [`build_exclude_covars`][fedcat.core.structures.build_exclude_covars].

        Args:
            chrom: Chromosome.

        Returns:
            Comma-separated covariates, or `"None"` if there are none.
        """
        return build_exclude_covars(
            chrom=chrom,
            pgs_prefix=self.pgs_prefix,
            chromosomes_no_pgs=self.chromosomes_no_pgs,
            exclude_covars=self.exclude_covars,
            exclude_covars_dict=self.exclude_covars_dict,
        )

    def step0_subset_dirname(self, cond_type: CondType) -> str:
        """Return the folder of a conditioning type's step-0 subsets, under the meta output folder.

        Named by the configured retest threshold, so subsets built at different
        thresholds never mix.

        Args:
            cond_type: Conditioning type.

        Returns:
            Relative folder, e.g. `step0_subsets/singles_p1e-06`.
        """
        p = self.cond_p_singles_retest if cond_type is CondType.SINGLEVAR else self.cond_p_agg_retest
        return f"step0_subsets/{cond_type.value}_p{p:g}"

    @property
    def step1_exclude_covars(self) -> str:
        """`--covarExcludeList` for step 1: `exclude_covars` only, or `"None"` if empty.

        The per-chromosome exclusions only apply to step 2, as step 1 fits every
        chromosome at once.
        """
        return build_exclude_covars(
            chrom="", pgs_prefix=None, chromosomes_no_pgs=(), exclude_covars=self.exclude_covars
        )

    def output_location(self, biobank: BioBank) -> str:
        """Return a biobank's results directory for this run.

        See [`compute_output_location`][fedcat.core.structures.compute_output_location].

        Args:
            biobank: Biobank.

        Returns:
            The results directory.
        """
        return compute_output_location(
            biobank=biobank,
            root_results_directory=self.cohorts[biobank].root_results_directory,
            run_prefix=self.run_prefix,
            trait_transform=self.trait_transform,
            trait_type=self.trait_type,
            aou_bucket=self.aou_bucket,
        )

    def dependencies_path(self, biobank: BioBank) -> str:
        """Return the `gs://` path of a biobank's staged workflow dependencies.

        Built from `aou_bucket` and the cohort's `dependencies_version` (the default
        version when there is no AoU cohort, e.g. for an AoU-hosted meta-analysis).

        Args:
            biobank: Biobank; only AoU is supported.

        Returns:
            The dependencies directory path.

        Raises:
            NotImplementedError: If `biobank` is not AoU.
            ValueError: If `aou_bucket` is unset.
        """
        if biobank is not BioBank.AOU:
            raise NotImplementedError(
                f"Workflow dependency staging is only implemented for AoU, not {biobank.value!r}"
            )
        if not self.aou_bucket:
            raise ValueError("aou_bucket is required to resolve the workflow dependencies path")
        cohort = self.cohorts.get(biobank)
        version = cohort.dependencies_version if cohort else DEFAULT_DEPENDENCIES_VERSION
        return compute_dependencies_path(self.aou_bucket, version)

    def validate(self) -> None:
        """Check field values, and fill `cohort_sizes` for every stratum.

        Raises:
            ValueError: Listing every invalid field.
        """
        errors: list[str] = []
        if not self.pheno_names:
            errors.append("pheno_names cannot be empty")
        bad_names = [p for p in self.pheno_names if not PHENO_NAME_RE.fullmatch(p)]
        if bad_names:
            errors.append(f"pheno_names may only contain letters, digits, '_', '.' and '-': {bad_names}")
        if len(set(self.pheno_names)) != len(self.pheno_names):
            errors.append("pheno_names contains duplicates")
        for label in ("step1_pheno_group_size", "step2_pheno_group_size"):
            if getattr(self, label) < 0:
                errors.append(f"{label} cannot be negative (0 = all phenotypes together)")
        for label in ("max_step2_singles_block_length", "max_step2_aggregates_block_length"):
            value = getattr(self, label)
            if value is not None and value < 100_000:
                errors.append(f"{label} must be at least 100000 (the variant block file resolution)")
        if self.regenie_version not in REGENIE_VERSIONS:
            errors.append(f"regenie_version {self.regenie_version!r} must be one of {list(REGENIE_VERSIONS)}")
        if not self.aggregate_freqs:
            errors.append("aggregate_freqs cannot be empty")
        else:
            for f in self.aggregate_freqs:
                try:
                    val = float(f)
                except (TypeError, ValueError):
                    errors.append(f"aggregate_freqs entry {f!r} is not a valid number")
                    continue
                # Must read exactly as regenie labels the bin in aggregate IDs (C++ default
                # formatting, i.e. Python's "g"), and as "0.x" so the bin suffix can be stripped
                if not (0.0 < val < 1.0) or f != format(val, "g") or not f.startswith("0."):
                    errors.append(
                        f"aggregate_freqs entry {f!r} must be in [0.0001, 1) and written as regenie "
                        f"labels it (e.g. '0.001', not '0.0010' or '1e-3')"
                    )
        if (strat := self.stratify_by_covar) is not None:
            if not strat.covar_name or any(ch.isspace() or ch == "," for ch in strat.covar_name):
                errors.append(
                    f"stratify_by_covar.covar_name {strat.covar_name!r} must be a single covariate column name"
                )
            if not strat.covar_values_names:
                errors.append("stratify_by_covar.covar_values_names cannot be empty")
            if "" in strat.covar_values_names:
                errors.append("stratify_by_covar.covar_values_names cannot contain an empty value")
            names = list(strat.covar_values_names.values())
            bad_names = [n for n in names if not PHENO_NAME_RE.fullmatch(n) or n == COVAR_STRATUM_COMBINED_NAME]
            if bad_names:
                errors.append(
                    "stratify_by_covar.covar_values_names names may only contain letters, digits, '_', '.' and '-', "
                    f"and cannot be {COVAR_STRATUM_COMBINED_NAME!r}: {bad_names}"
                )
            if len(set(names)) != len(names):
                errors.append("stratify_by_covar.covar_values_names contains duplicate names")
        if self.interaction and not INTERACTION_RE.fullmatch(self.interaction):
            errors.append(
                f"interaction {self.interaction!r} must be a covariate column name, "
                "optionally followed by a baseline level in brackets, e.g. 'BMI' or 'smoking[0]'"
            )
        elif strat is not None and self.interaction.split("[")[0] == strat.covar_name:
            errors.append(
                f"interaction {self.interaction!r} cannot be the stratifying covariate; "
                "it is constant within each covariate stratum"
            )
        # Excluded covariate names are comma-joined for regenie's --covarExcludeList.
        for label, names in [
            ("exclude_covars", self.exclude_covars),
            *((f"exclude_covars_dict[{c!r}]", v) for c, v in self.exclude_covars_dict.items()),
        ]:
            for name in names:
                if "," in name or any(ch.isspace() for ch in name):
                    errors.append(f"{label} entry {name!r} cannot contain ',' or whitespace")
        for biobank, cohort in self.cohorts.items():
            if not cohort.strata:
                errors.append(f"cohort {biobank.value} has no strata")
            bad_strata = [s for s in cohort.strata if not STRATUM_NAME_RE.fullmatch(s)]
            if bad_strata:
                errors.append(
                    f"cohort {biobank.value} strata may only contain letters, digits and '-': {bad_strata}"
                )
            if len(set(cohort.strata)) != len(cohort.strata):
                errors.append(f"cohort {biobank.value} strata contains duplicates")
            # Passed to the workflow scripts as an unquoted positional argument
            if not cohort.cat_vars or any(ch.isspace() for ch in cohort.cat_vars):
                errors.append(
                    f"cohort {biobank.value} cat_vars must be a comma-separated list or NONE, "
                    f"without spaces: {cohort.cat_vars!r}"
                )
            for label in ("max_tasks_per_workflow", "max_concurrent_jobs"):
                if getattr(cohort, label) < 0:
                    errors.append(f"cohort {biobank.value} {label} cannot be negative (0 = no limit)")
            for kind, resources in cohort.workflows.resources.items():
                where = f"cohort {biobank.value} resources.{kind.value}"
                for label in ("max_retries", "preemptible"):
                    value = getattr(resources, label)
                    if value is not None and value < 0:
                        errors.append(f"{where} {label} cannot be negative")
                if resources.priority is not None and resources.priority not in DX_PRIORITIES:
                    errors.append(f"{where} priority must be one of {sorted(DX_PRIORITIES)}")
            missing_sizes = [s for s in cohort.strata if s not in cohort.cohort_sizes]
            if missing_sizes:
                logger.warning(
                    f"Biobank {biobank.value} missing sizes for strata: {missing_sizes} (default=100,000 used)"
                )
            unknown_sizes = [s for s in cohort.cohort_sizes if s not in cohort.strata]
            if unknown_sizes:
                logger.warning(f"Biobank {biobank.value} has sizes for strata it does not run: {unknown_sizes}")
            cohort.cohort_sizes = {s: cohort.cohort_sizes.get(s, 100_000) for s in cohort.strata}
        if errors:
            raise ValueError("Configuration validation failed: " + "; ".join(errors))


# ── Runtime-only structures ───────────────────────────────────────────────────


@define
class Phenotype:
    """One phenotype of a [`DiscoveryConfig`][fedcat.core.structures.DiscoveryConfig]."""

    name: str
    trait_type: TraitType
    transform: Transform

    @classmethod
    def from_config(cls, config: DiscoveryConfig, pheno_name: str) -> Phenotype:
        """Create a phenotype with the config's trait type and transform.

        Args:
            config: Run config.
            pheno_name: Phenotype name.

        Returns:
            The phenotype.
        """
        return cls(
            name=pheno_name,
            trait_type=config.trait_type,
            transform=config.trait_transform,
        )


def _drop_covar(cat_vars: str, covar: str) -> str:
    """Remove `covar` from a comma-separated `cat_vars`, returning `"NONE"` if none are left."""
    kept = [name for name in (n.strip() for n in cat_vars.split(",")) if name and name != covar]
    if not kept:
        return "NONE"
    return ",".join(kept)


def build_exclude_covars(
    chrom: str,
    pgs_prefix: str | None,
    chromosomes_no_pgs: Sequence[str],
    exclude_covars: Sequence[str] = (),
    exclude_covars_dict: Mapping[str, Sequence[str]] | None = None,
) -> str:
    """Build REGENIE's `--covarExcludeList` for one chromosome.

    Combines, in order and without duplicates: the PGS column, `exclude_covars`, and
    `exclude_covars_dict[chrom]`.

    Args:
        chrom: Chromosome.
        pgs_prefix: PGS covariate prefix; the column `{pgs_prefix}{chrom}` is excluded
            unless `chrom` is in `chromosomes_no_pgs`. None for no PGS covariate.
        chromosomes_no_pgs: Chromosomes without a PGS covariate.
        exclude_covars: Covariates excluded on every chromosome.
        exclude_covars_dict: Further covariates excluded per chromosome.

    Returns:
        Comma-separated covariates, or `"None"`, which the workflows take to mean no
        `--covarExcludeList`.
    """
    names: list[str] = []
    if pgs_prefix and chrom not in chromosomes_no_pgs:
        names.append(f"{pgs_prefix}{chrom}")
    names.extend(exclude_covars)
    names.extend((exclude_covars_dict or {}).get(chrom, ()))
    deduped = list(dict.fromkeys(n for n in names if n))
    return ",".join(deduped) if deduped else "None"


def compute_output_location(
    biobank: BioBank,
    root_results_directory: str,
    run_prefix: str,
    trait_transform: Transform,
    trait_type: TraitType,
    aou_bucket: str | None = None,
) -> str:
    """Return a cohort's results directory for a run.

    Args:
        biobank: Biobank.
        root_results_directory: Root results directory in the cohort config.
        run_prefix: Run prefix.
        trait_transform: Phenotype transform.
        trait_type: Trait type; names the directory of `RAW` runs.
        aou_bucket: Workspace bucket, required for AoU.

    Returns:
        `<bucket>/<root>/<run_prefix>_<transform>_wgs` for AoU, or
        `/<root>/<run_prefix>_<transform>_wgs` for UKB, where `<transform>` is `rint`,
        `mcc`, `raw` (continuous) or `binary`.

    Raises:
        ValueError: If `biobank` is AoU and `aou_bucket` is unset.
    """
    if trait_transform == Transform.RINT:
        outputdir = f"{run_prefix}_rint_wgs"
    elif trait_transform == Transform.MCC:
        outputdir = f"{run_prefix}_mcc_wgs"
    elif trait_type == TraitType.CONTINUOUS:
        outputdir = f"{run_prefix}_raw_wgs"
    else:
        outputdir = f"{run_prefix}_binary_wgs"

    if biobank == BioBank.AOU:
        if not aou_bucket:
            raise ValueError("aou_bucket is required for AOU cohorts")
        return f"{aou_bucket.removeprefix('gs://')}/{root_results_directory.lstrip('/')}/{outputdir}"
    else:  # UKB
        return f"/{root_results_directory.lstrip('/')}/{outputdir}"


def stratum_descrip(biobank: BioBank, stratum: str) -> str:
    """Return the `{stratum}_{biobank}` label of a cohort stratum, e.g. `eur_aou`."""
    return f"{stratum}_{biobank.value}"


# Config fields whose default would change the results of a run already in progress.
# A persisted config that predates one gets the value that reproduces its earlier
# steps; resuming with the new default is then a config difference to resolve.
LEGACY_CONFIG_VALUES: dict[str, Any] = {"acat_singles_stratum_rare": False}


def pin_legacy_config(raw_config: dict[str, Any]) -> list[str]:
    """Fill fields missing from a persisted config with their legacy values.

    Args:
        raw_config: Persisted config as parsed JSON; updated in place.

    Returns:
        Names of the fields filled.
    """
    missing = [k for k in LEGACY_CONFIG_VALUES if k not in raw_config]
    for k in missing:
        raw_config[k] = LEGACY_CONFIG_VALUES[k]
    return missing


def make_conv():
    """Create the cattrs converter for config and state JSON, with strict scalar types.

    A value of the wrong type raises ValueError, as it is a config error.

    Returns:
        The converter.
    """
    conv = make_converter()

    @conv.register_structure_hook
    def string_hook(value, _cl) -> str:
        """Accept only a string."""
        if not isinstance(value, str):
            raise ValueError(f"{value!r} in config file should be a string")  # noqa: TRY004
        return value

    @conv.register_structure_hook
    def list_hook(value, _cl) -> list:
        """Accept only a list."""
        if not isinstance(value, list):
            raise ValueError(f"{value!r} in config file should be a list")  # noqa: TRY004
        return value

    @conv.register_structure_hook
    def float_hook(value, _cl) -> float:
        """Accept only a float, so `1` is rejected where `1.0` is meant."""
        if not isinstance(value, float):
            raise ValueError(f"{value!r} in config file should be a float")  # noqa: TRY004
        return value

    def bool_hook(value, _cl) -> bool:
        """Accept only a boolean."""
        if not isinstance(value, bool):
            raise ValueError(f"{value!r} in config file should be a boolean")  # noqa: TRY004
        return value

    def int_hook(value, _cl) -> int:
        """Accept only an integer."""
        if not isinstance(value, int):
            raise ValueError(f"{value!r} in config file should be an integer")  # noqa: TRY004
        return value

    # Exact-type matches: a class-based int hook would also catch IntEnum subclasses
    # (JobState, StepState) and structure them as plain ints
    conv.register_structure_hook_func(lambda t: t is bool, bool_hook)
    conv.register_structure_hook_func(lambda t: t is int, int_hook)

    return conv
