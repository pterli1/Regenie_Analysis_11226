"""Preflight checks of a run's config, authentication, workflows and input files.

A check is an `async def check_*(ctx: CheckContext) -> list[Finding]` that returns its
findings rather than raising or logging. Run them with
[`run_checks`][fedcat.core.checks.run_checks] and report them with
[`log_findings`][fedcat.core.checks.log_findings].

Checks run by scope, in the order `config`, `auth`, `workflows`, `remote_files`; an
error in `config` or `auth` skips the later scopes. `--skip-validation` skips
`workflows` and `remote_files`.
"""

from __future__ import annotations

import asyncio
from enum import Enum
from pathlib import Path
from typing import TYPE_CHECKING

from attrs import define, field
from loguru import logger

from fedcat.core.structures import (
    BIOBANK_PROFILES,
    BioBank,
    EgressPolicy,
    MetaMode,
    Phase,
    RunEnv,
)
from fedcat.platforms.base import PathSyntaxError, PlatformError

if TYPE_CHECKING:
    from fedcat.core.structures import DiscoveryConfig
    from fedcat.platforms.base import Platform


# ── types ─────────────────────────────────────────────────────────────────────


class Severity(Enum):
    """Severity of a finding; any error stops the run."""

    ERROR = "error"
    WARNING = "warning"


# Check scopes, in run order
CONFIG = "config"
AUTH = "auth"
WORKFLOWS = "workflows"
REMOTE_FILES = "remote_files"


@define(frozen=True)
class Finding:
    """A problem found by a check."""

    severity: Severity
    """Severity."""
    scope: str
    """Scope of the check that found it."""
    label: str
    """What it concerns, e.g. a config field."""
    message: str
    """Description of the problem."""


@define(frozen=True)
class RunFlags:
    """Command-line flags that change which checks run."""

    dry_run: bool = False
    """`--dry-run`: platforms connect without their job services where possible."""
    skip_validation: bool = False
    """`--skip-validation`: skip the `workflows` and `remote_files` scopes."""


@define(frozen=True)
class CheckContext:
    """What the checks examine."""

    config: DiscoveryConfig
    """Run config."""
    platforms: dict[BioBank, Platform]
    """Platform per biobank."""
    phases_to_run: list[Phase]
    """Phases to be run."""
    variant_density_blockfile: Path | None = None
    """Variant density block file."""
    meta_mode: MetaMode = MetaMode.ALL
    """Meta-analysis mode."""
    flags: RunFlags = field(factory=RunFlags)
    """Command-line flags."""


def _needs_aggregates(ctx: CheckContext) -> bool:
    """Return whether a phase using the cohorts' aggregate annotations will run."""
    return Phase.STEP2_AGG in ctx.phases_to_run or Phase.COND_AGG in ctx.phases_to_run


def _runs_meta_jobs(ctx: CheckContext) -> bool:
    """Return whether a phase running `meta_analysis_config` jobs (meta-analysis or step-0 subset) will run."""
    return Phase.META in ctx.phases_to_run or Phase.STEP0_SUBSET in ctx.phases_to_run


def _needs_meta_aggregates(ctx: CheckContext) -> bool:
    """Return whether the meta-analysis or step-0 subset will handle burden/ACAT results, which need annotations."""
    return _runs_meta_jobs(ctx) and ctx.meta_mode.runs_aggregates


# ── config scope (pure, no I/O) ───────────────────────────────────────────────


async def check_egress_policy(ctx: CheckContext) -> list[Finding]:
    """Check fedcat runs inside the environment of every cohort that forbids data egress."""
    findings: list[Finding] = []
    for bb in ctx.config.cohorts:
        if BIOBANK_PROFILES[bb].egress is EgressPolicy.FORBIDDEN and ctx.config.local_env is not RunEnv.AOU:
            findings.append(
                Finding(
                    Severity.ERROR,
                    CONFIG,
                    f"{bb.value}.egress",
                    f"cohort {bb.value!r} forbids data egress, so fedcat must run inside its "
                    f"environment with local_env='aou' (got {ctx.config.local_env.value!r}). "
                    "Downloading its results to another platform is not permitted.",
                )
            )
    return findings


async def check_meta_config_present(ctx: CheckContext) -> list[Finding]:
    """Check `meta_analysis_config`, and the workflow of each of its phases being run, are set."""
    if not _runs_meta_jobs(ctx):
        return []
    meta_cfg = ctx.config.meta_analysis_config
    if meta_cfg is None:
        return [
            Finding(
                Severity.ERROR,
                CONFIG,
                "meta_analysis_config",
                "meta_analysis_config must be set in the config to run the meta_analysis / step0_subset phases",
            )
        ]
    return [
        Finding(
            Severity.ERROR,
            CONFIG,
            f"meta_analysis_config.{field_name}",
            f"meta_analysis_config.{field_name} must be set to run the {phase.value} phase",
        )
        for phase, field_name, value in (
            (Phase.META, "meta_workflow", meta_cfg.meta_workflow),
            (Phase.STEP0_SUBSET, "subset_workflow", meta_cfg.subset_workflow),
        )
        if phase in ctx.phases_to_run and not value
    ]


async def check_step0_between_meta_and_cond(ctx: CheckContext) -> list[Finding]:
    """Check `step0_subset` runs between a meta-analysis and conditioning run together.

    Otherwise conditioning would start from a subset older than the meta-analysis.
    """
    phases = ctx.phases_to_run
    cond = [p for p in (Phase.COND_SINGLES, Phase.COND_AGG) if p in phases]
    if Phase.META in phases and cond and Phase.STEP0_SUBSET not in phases:
        return [
            Finding(
                Severity.ERROR,
                CONFIG,
                "phases_to_run",
                f"{Phase.META.value} and {', '.join(p.value for p in cond)} need "
                f"{Phase.STEP0_SUBSET.value} between them in --phases-to-run",
            )
        ]
    return []


async def check_meta_mode(ctx: CheckContext) -> list[Finding]:
    """Check the meta mode can run and covers the results the conditioning phases being run need."""
    if not _runs_meta_jobs(ctx):
        return []
    mode = ctx.meta_mode
    have_cond_file = ctx.config.has_cond_file
    if mode is MetaMode.COND_SINGLES_ONLY and not have_cond_file:
        return [
            Finding(
                Severity.ERROR,
                CONFIG,
                "meta_mode",
                f"meta mode {mode.value!r} needs a cond_file in every cohort",
            )
        ]
    stages = mode.stages(have_cond_file)
    findings: list[Finding] = []
    for phase, ok, needs in (
        (Phase.COND_SINGLES, stages.uncond_singles, "unconditional singles"),
        (Phase.COND_AGG, stages.aggregates, "aggregates"),
    ):
        if phase in ctx.phases_to_run and not ok:
            findings.append(
                Finding(
                    Severity.ERROR,
                    CONFIG,
                    "meta_mode",
                    f"meta mode {mode.value!r} does not cover the {needs} "
                    f"that {phase.value} needs; use --meta-mode auto or drop {phase.value}",
                )
            )
    return findings


async def check_cond_suffix(ctx: CheckContext) -> list[Finding]:
    """Check `cond_suffix` is set when a conditioning phase runs."""
    needs_cond = Phase.COND_SINGLES in ctx.phases_to_run or Phase.COND_AGG in ctx.phases_to_run
    if needs_cond and not ctx.config.cond_suffix:
        return [
            Finding(
                Severity.ERROR,
                CONFIG,
                "cond_suffix",
                "cond_suffix must be set to run the conditioning phases",
            )
        ]
    return []


async def check_covar_stratification(ctx: CheckContext) -> list[Finding]:
    """Check a covariate stratum's config excludes the stratifying covariate from the covariates.

    It is constant within the stratum, so must be in `exclude_covars` and not in any
    `cat_vars`, as [`covar_stratum_configs`][fedcat.core.structures.DiscoveryConfig.covar_stratum_configs]
    sets it.
    """
    cfg = ctx.config
    strat = cfg.stratify_by_covar
    if strat is None or not cfg.covar_stratum:
        return []

    findings: list[Finding] = []
    if strat.covar_name not in cfg.exclude_covars:
        findings.append(
            Finding(
                Severity.ERROR,
                CONFIG,
                "stratify_by_covar.covar_name",
                f"stratifying covariate {strat.covar_name!r} is not in exclude_covars; it is constant "
                f"within the {cfg.covar_stratum!r} covariate stratum, so it must be excluded from the model",
            )
        )
    findings += [
        Finding(
            Severity.ERROR,
            CONFIG,
            f"{bb.value}.cat_vars",
            f"stratifying covariate {strat.covar_name!r} is still listed in cohort {bb.value!r}'s "
            "cat_vars; it cannot be both a categorical covariate and excluded",
        )
        for bb, cohort in cfg.cohorts.items()
        if strat.covar_name in [n.strip() for n in cohort.cat_vars.split(",")]
    ]
    return findings


async def check_cohort_annotations(ctx: CheckContext) -> list[Finding]:
    """Check every cohort has annotations when an aggregate phase runs."""
    if not _needs_aggregates(ctx):
        return []
    return [
        Finding(
            Severity.ERROR,
            CONFIG,
            f"{bb.value}.regenie_annotations",
            f"cohort {bb.value!r} has no regenie_annotations, which are required to run "
            "step2_aggregates / conditioning_aggregates",
        )
        for bb, cohort in ctx.config.cohorts.items()
        if not cohort.regenie_annotations
    ]


async def check_meta_annotations(ctx: CheckContext) -> list[Finding]:
    """Check the meta-analysis has at most one annotation source, and one when it runs aggregates."""
    meta_cfg = ctx.config.meta_analysis_config
    if meta_cfg is None:
        return []
    has_schema_dir = bool(meta_cfg.schema_directory)
    has_annots = bool(meta_cfg.regenie_annotations)
    if has_schema_dir and has_annots:
        return [
            Finding(
                Severity.ERROR,
                CONFIG,
                "meta_analysis_config.annotations",
                "cannot set both meta_analysis_config.schema_directory and "
                "meta_analysis_config.regenie_annotations",
            )
        ]
    if not has_schema_dir and not has_annots and _needs_meta_aggregates(ctx):
        return [
            Finding(
                Severity.ERROR,
                CONFIG,
                "meta_analysis_config.annotations",
                "must set one of meta_analysis_config.schema_directory or "
                f"meta_analysis_config.regenie_annotations to meta-analyse aggregates "
                f"(meta mode {ctx.meta_mode.value!r})",
            )
        ]
    return []


# ── auth scope ────────────────────────────────────────────────────────────────


async def check_auth(ctx: CheckContext) -> list[Finding]:
    """Check each cohort's platform connects, with its job service unless this is a dry run."""
    findings: list[Finding] = []
    need_submit = not ctx.flags.dry_run
    for bb, platform in ctx.platforms.items():
        try:
            await platform.connect(need_submit=need_submit)
        except PlatformError as e:
            findings.append(Finding(Severity.ERROR, AUTH, f"{bb.value}.auth", str(e)))
    return findings


# ── workflows scope ───────────────────────────────────────────────────────────


def _required_targets(ctx: CheckContext) -> dict[BioBank, list[str]]:
    """Return the workflows or applets the phases need, per biobank they run on."""
    phases = set(ctx.phases_to_run)
    out: dict[BioBank, list[str]] = {}
    for bb, cohort in ctx.config.cohorts.items():
        wf = cohort.workflows
        phase_map = {
            Phase.STEP1: wf.regenie_step1,
            Phase.STEP2_SINGLES: wf.regenie_step2_singlevar,
            Phase.STEP2_AGG: wf.regenie_step2_aggregate,
            Phase.COND_SINGLES: wf.conditioning_singlevar,
            Phase.COND_AGG: wf.conditioning_aggregate,
        }
        targets = [w for ph, w in phase_map.items() if ph in phases]
        if targets:
            out.setdefault(bb, []).extend(targets)
    meta_cfg = ctx.config.meta_analysis_config
    if meta_cfg is not None:
        meta_bb = BioBank.AOU if meta_cfg.location is RunEnv.AOU else BioBank.UKB
        for phase, target in (
            (Phase.META, meta_cfg.meta_workflow),
            (Phase.STEP0_SUBSET, meta_cfg.subset_workflow),
        ):
            if phase in phases and target and meta_bb in ctx.platforms:
                out.setdefault(meta_bb, []).append(target)
    return out


async def check_workflows(ctx: CheckContext) -> list[Finding]:
    """Check the workflows or applets the phases need exist on their platforms."""
    findings: list[Finding] = []
    for bb, targets in _required_targets(ctx).items():
        problems = await ctx.platforms[bb].check_targets(targets)
        findings.extend(Finding(Severity.ERROR, WORKFLOWS, f"{bb.value}.workflow", p) for p in problems)
    return findings


async def check_blockfile(ctx: CheckContext) -> list[Finding]:
    """Check the variant density block file exists when a step 2 phase runs."""
    needs_step2 = Phase.STEP2_SINGLES in ctx.phases_to_run or Phase.STEP2_AGG in ctx.phases_to_run
    bf = ctx.variant_density_blockfile
    if needs_step2 and bf is not None and not bf.exists():
        return [
            Finding(
                Severity.ERROR,
                WORKFLOWS,
                "variant_density_blockfile",
                f"required for step2 but not found: {bf}",
            )
        ]
    return []


# ── remote_files scope ────────────────────────────────────────────────────────


def _build_remote_checks(
    ctx: CheckContext,
) -> tuple[set[tuple[str, str, BioBank]], list[tuple[str, BioBank, list[list[str]]]]]:
    """Return the remote files the run reads, to check they exist.

    Returns:
        `(checks, alt_checks)`: `checks` is `{(label, path, biobank)}` of files that must
        all exist; `alt_checks` is `[(label, biobank, alternatives)]`, satisfied when
        every file of one alternative exists (PLINK1 or PLINK2 step 1 input).
    """
    from fedcat.stats.subset_results import _PARQUET_SCHEMA_FILES

    config = ctx.config
    checks: set[tuple[str, str, BioBank]] = set()
    alt_checks: list[tuple[str, BioBank, list[list[str]]]] = []
    needs_aggregates = _needs_aggregates(ctx)

    for bb, cohort in config.cohorts.items():
        lbl = bb.value
        for name, path in [("pheno_file", cohort.pheno_file), ("cov_file", cohort.cov_file)]:
            if path:
                checks.add((f"{lbl}.{name}", path, bb))
        if cohort.cond_file:
            checks.add((f"{lbl}.cond_file", cohort.cond_file, bb))

        for chrom in config.chromosomes_to_run:
            resolved = cohort.pfile_template.replace("{chr}", chrom)
            for ext in (".pgen", ".psam", ".pvar"):
                checks.add((f"{lbl}.pfile chr{chrom}{ext}", f"{resolved}{ext}", bb))

        for stratum in cohort.strata:
            resolved = cohort.step1_input_template.replace("{stratum}", stratum)
            alt_checks.append(
                (
                    f"{lbl}.step1_input {stratum}",
                    bb,
                    [
                        [f"{resolved}{ext}" for ext in exts]
                        for exts in ((".bed", ".bim", ".fam"), (".pgen", ".pvar", ".psam"))
                    ],
                )
            )

        if needs_aggregates:
            for annot in cohort.regenie_annotations:
                for chrom in config.chromosomes_to_run:
                    for fname, fval in [
                        ("mask_file", annot.mask_file),
                        ("annot_file", annot.annot_file),
                        ("setlist_file", annot.setlist_file),
                    ]:
                        resolved = fval.replace("{chr}", chrom)
                        checks.add((f"{lbl}.annot[{annot.name}].{fname} chr{chrom}", resolved, bb))

        if BIOBANK_PROFILES[bb].platform_kind == "gcs_cromwell" and config.aou_bucket:
            from fedcat.platforms.gcs_cromwell.workflow_setup import dependency_filenames

            dep_dir = config.dependencies_path(bb)
            for fname in dependency_filenames(bb, cond_engine=config.cond_engine):
                checks.add((f"{lbl}.dependencies/{fname}", f"{dep_dir}/{fname}", bb))

    meta_cfg = config.meta_analysis_config
    if meta_cfg is not None and _needs_meta_aggregates(ctx):
        meta_bb = BioBank.AOU if meta_cfg.location is RunEnv.AOU else BioBank.UKB
        if meta_cfg.schema_directory:
            for chrom in config.chromosomes_to_run:
                resolved_dir = meta_cfg.schema_directory.replace("{chr}", chrom).rstrip("/")
                for schema_file in _PARQUET_SCHEMA_FILES:
                    checks.add(
                        (
                            f"meta.schema_directory chr{chrom}/{schema_file}",
                            f"{resolved_dir}/{schema_file}",
                            meta_bb,
                        )
                    )
        elif meta_cfg.regenie_annotations:
            for annot in meta_cfg.regenie_annotations:
                for chrom in config.chromosomes_to_run:
                    for fname, fval in [
                        ("mask_file", annot.mask_file),
                        ("annot_file", annot.annot_file),
                        ("setlist_file", annot.setlist_file),
                    ]:
                        resolved = fval.replace("{chr}", chrom)
                        checks.add((f"meta.annot[{annot.name}].{fname} chr{chrom}", resolved, meta_bb))

    return checks, alt_checks


async def check_remote_files(ctx: CheckContext) -> list[Finding]:
    """Check every remote input file the run reads exists, checking up to 16 at once."""
    checks, alt_checks = _build_remote_checks(ctx)
    sem = asyncio.Semaphore(16)

    async def _missing(path: str, bb: BioBank) -> str | None:
        """Return the problem with a path (missing or malformed), or None if it exists."""
        platform = ctx.platforms.get(bb)
        if platform is None:
            return f"no platform configured for {bb.value}"
        async with sem:
            try:
                return None if await platform.exists(path) else "not found"
            except PathSyntaxError as e:
                return f"malformed reference ({e})"
            except PlatformError as e:
                return str(e)

    findings: list[Finding] = []

    direct = list(checks)
    direct_results = await asyncio.gather(*(_missing(p, bb) for _, p, bb in direct))
    for (label, path, _bb), problem in zip(direct, direct_results):
        if problem is not None:
            findings.append(Finding(Severity.ERROR, REMOTE_FILES, label, f"{path}: {problem}"))

    # Alternative filesets
    alt_pairs = sorted({(p, bb) for _, bb, alts in alt_checks for alt in alts for p in alt})
    alt_results = await asyncio.gather(*(_missing(p, bb) for p, bb in alt_pairs))
    missing_alt = {pair for pair, problem in zip(alt_pairs, alt_results) if problem is not None}
    for label, bb, alts in alt_checks:
        if any(all((p, bb) not in missing_alt for p in alt) for alt in alts):
            continue
        wanted = " or ".join(", ".join(alt) for alt in alts)
        findings.append(
            Finding(Severity.ERROR, REMOTE_FILES, label, f"no complete fileset, need all of {wanted}")
        )

    return sorted(findings, key=lambda f: f.label)


# ── scope registration ────────────────────────────────────────────────────────

CONFIG_CHECKS = [
    check_egress_policy,
    check_meta_config_present,
    check_step0_between_meta_and_cond,
    check_meta_mode,
    check_cond_suffix,
    check_covar_stratification,
    check_cohort_annotations,
    check_meta_annotations,
]
AUTH_CHECKS = [check_auth]
WORKFLOW_CHECKS = [check_workflows, check_blockfile]
REMOTE_CHECKS = [check_remote_files]


# ── runner ────────────────────────────────────────────────────────────────────


async def _run_scope(scope: str, checks, ctx: CheckContext) -> list[Finding]:
    """Run a scope's checks, turning a check that raises into an error finding."""
    findings: list[Finding] = []
    for fn in checks:
        try:
            findings.extend(await fn(ctx))
        except Exception as exc:  # noqa: BLE001 — a bug in a check must surface, never pass silently
            findings.append(Finding(Severity.ERROR, scope, fn.__name__, f"check crashed: {exc!r}"))
    return findings


def _has_error(findings: list[Finding]) -> bool:
    """Return whether any finding is an error."""
    return any(f.severity is Severity.ERROR for f in findings)


async def run_checks(ctx: CheckContext) -> list[Finding]:
    """Run the checks scope by scope, stopping after a `config` or `auth` error.

    Args:
        ctx: What to check.

    Returns:
        All findings.
    """
    findings = await _run_scope(CONFIG, CONFIG_CHECKS, ctx)
    if _has_error(findings):
        return findings

    auth = await _run_scope(AUTH, AUTH_CHECKS, ctx)
    findings += auth
    if _has_error(auth):
        return findings

    if not ctx.flags.skip_validation:
        findings += await _run_scope(WORKFLOWS, WORKFLOW_CHECKS, ctx)
        findings += await _run_scope(REMOTE_FILES, REMOTE_CHECKS, ctx)

    return findings


def log_findings(findings: list[Finding]) -> bool:
    """Log each finding at its severity.

    Args:
        findings: Findings from [`run_checks`][fedcat.core.checks.run_checks].

    Returns:
        True if any finding is an error.
    """
    errors = [f for f in findings if f.severity is Severity.ERROR]
    warnings = [f for f in findings if f.severity is Severity.WARNING]
    for f in warnings:
        logger.warning(f"[{f.scope}] {f.label}: {f.message}")
    if errors:
        logger.error(f"Preflight checks failed — {len(errors)} error(s):")
        for f in errors:
            logger.error(f"  [{f.scope}] {f.label}: {f.message}")
    return bool(errors)
