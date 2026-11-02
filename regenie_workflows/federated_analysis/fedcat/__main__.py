"""The `fedcat` command line: argument parsing and the subcommands.

The options of each subcommand are documented in their argparse help, from which
`tools/docs/gen_docs.py` writes `docs/options_reference.md` and the documentation's options pages.
"""

from __future__ import annotations

# Sets the BLAS thread variables, so must be imported before anything that loads numpy
import fedcat._blas_threads  # noqa: F401

# isort: split
import argparse
import asyncio
import contextlib
import json
import math
import sys
from datetime import datetime, timedelta
from pathlib import Path

from loguru import logger

from fedcat import __version__
from fedcat.core.engine import STATE_VERSION, DiscoveryRun
from fedcat.core.logging import STAGE_VAR, STATE_WRITE_QUEUE, setup_run_logging
from fedcat.core.structures import (
    COVAR_STRATUM_COMBINED,
    COVAR_STRATUM_COMBINED_NAME,
    DEFAULT_DEPENDENCIES_VERSION,
    PHASE_ORDERED,
    PHASE_STRINGS_ORDERED,
    BioBank,
    DiscoveryConfig,
    JobState,
    MetaMode,
    Phase,
    RunEnv,
    WorkflowBackend,
    make_conv,
)


def _load_discovery_state(state_path: Path, conv) -> DiscoveryRun:
    """Load a phase's state file, exiting if it was saved by an incompatible fedcat version.

    The version is read from the raw JSON, as structuring would fill in the current
    version for a state file without one.
    """
    raw = json.loads(state_path.read_text())
    persisted_version = raw.get("state_version", "UNVERSIONED")
    if persisted_version != STATE_VERSION:
        logger.error(
            f"State file {state_path} has schema version {persisted_version}, but this fedcat"
            f"(v{__version__}) expects schemaversion {STATE_VERSION}. The persisted state format changed: "
            f"start a fresh run (remove the old state file); manually edit state file; "
            f"or downgrade fedcat to the version that created the state file. "
            f"Compatible versions:"
            f"\n\tschema version 5\tfedcat >= 0.8.0"
            f"\n\tschema version 4\tfedcat 0.7.0 - 0.7.2"
            f"\n\tschema version 3\tfedcat 0.5.0 - 0.6.3"
            f"\n\tschema version 2\tfedcat 0.4.0 - 0.4.1"
            f"\n\tUNVERSIONED\tfedcat <= 0.3.1 (try to match based on file age)"
        )
        sys.exit(1)
    _pin_persisted_config(raw.get("config", {}), state_path)
    return conv.structure(raw, DiscoveryRun)


def _pin_persisted_config(raw_config: dict, state_path: Path) -> None:
    """Give a persisted config the legacy value of options added since, warning for each."""
    from fedcat.core.structures import pin_legacy_config

    for key in pin_legacy_config(raw_config):
        logger.warning(
            f"{state_path} predates {key!r}; keeping its earlier behaviour ({key}=false). "
            f"Set it to false in the config to resume, or re-run from the meta-analysis to use true."
        )


def _set_aside_later_phases(run_dir: Path, phase: Phase) -> None:
    """Rename the state files of phases after `phase` to `.old`, so a resume does not start from them."""
    for later_phase in PHASE_ORDERED[PHASE_ORDERED.index(phase) + 1 :]:
        state = run_dir / f"discovery_state_{later_phase.value}.json"
        if state.exists():
            state.replace(state.with_name(f"{state.name}.old"))
            logger.info(f"Moved {state.name} to {state.name}.old")


def _persisted_workspace_bucket(run_dir: Path, existing_phases: list[Phase]) -> tuple[str, str] | None:
    """Return the `(aou_bucket, aou_bucket_resource_id)` saved in the latest phase state that has one.

    The bucket chosen on a run's first invocation is reused rather than asked for again.
    """
    for phase in reversed(existing_phases):
        state_path = run_dir / f"discovery_state_{phase.value}.json"
        try:
            cfg = json.loads(state_path.read_text()).get("config", {})
        except (OSError, json.JSONDecodeError):
            continue
        if bucket := cfg.get("aou_bucket"):
            return bucket, cfg.get("aou_bucket_resource_id", "")
    return None


async def _abort_all(config_obj: DiscoveryConfig, run_dir: Path, conv) -> None:
    """`--abort-all`: abort every running job recorded in the run's state files, after confirmation."""
    from fedcat.core.abort import abort_jobs, running_jobs
    from fedcat.core.tracking import run_units
    from fedcat.platforms import build_platforms

    platforms = build_platforms(config_obj)
    jobs = running_jobs([u.state_file for u in run_units(config_obj, run_dir, platforms)], conv)
    if not jobs:
        logger.info(f"No running jobs recorded in the state for {config_obj.run_prefix!r}; nothing to abort.")
        return
    logger.info(f"{len(jobs)} running job(s) recorded for {config_obj.run_prefix!r}:")
    for job in jobs:
        logger.info(f"  {job.biobank.value:<4} {job.description} ({job.id})")
    await logger.complete()
    if input(f"Abort all {len(jobs)} job(s)? y/N: ").strip().lower() != "y":
        logger.info("Abort cancelled.")
        return
    n_aborted = await abort_jobs(jobs, platforms)
    logger.info(f"Abort requested for {n_aborted}/{len(jobs)} job(s).")


async def _report(config_obj: DiscoveryConfig, run_dir: Path, conv) -> None:
    """`--report`: log job states, task counts and costs of every phase and conditioning step."""
    from fedcat.core.tracking import (
        connect_platforms,
        format_job_states,
        format_task_counts,
        load_jobs,
        run_units,
        task_counts_by_biobank,
        unit_costs,
    )
    from fedcat.platforms import build_platforms

    platforms = build_platforms(config_obj)
    units = [(unit, load_jobs(unit.state_file, conv)) for unit in run_units(config_obj, run_dir, platforms)]
    if not units:
        logger.info(f"No state recorded for {config_obj.run_prefix!r}; nothing to report.")
        return
    submitted = {j.biobank for _, jobs in units for j in jobs if j.id}
    connected = await connect_platforms(platforms, submitted, storage_is_enough=True)

    lines = [f"Report for {config_obj.run_prefix!r}:"]
    totals: dict[str, float] = {}
    estimated = False
    for unit, jobs in units:
        lines.append(f"  {unit.label}: {format_job_states(jobs)}")
        task_counts = await task_counts_by_biobank(jobs, connected)
        costs = await unit_costs(unit, jobs, connected)
        for bb in sorted(task_counts.keys() | costs.keys(), key=lambda b: b.value):
            expected = sum(j.expected_tasks for j in jobs if j.biobank is bb and j.id)
            tasks = format_task_counts(task_counts[bb], expected) if bb in task_counts else "unavailable"
            cost = costs.get(bb)
            if cost is None:
                cost_str = "unavailable"
            else:
                cost_str = f"{cost.amount:,.2f} {cost.currency} ({cost.basis})"
                totals[cost.currency] = totals.get(cost.currency, 0.0) + cost.amount
                estimated |= cost.estimated
            lines.append(f"    {bb.value:<4} {tasks}; cost: {cost_str}")
    if totals:
        lines.append(
            "  total cost: " + ", ".join(f"{amount:,.2f} {cur}" for cur, amount in sorted(totals.items()))
        )
    if estimated:
        lines.append(
            "  Estimated costs price compute and disk at list rates rather than reporting what was billed; "
            "each line says what its figure was derived from."
        )
    logger.info("\n".join(lines))


def _parse_when(tokens: list[str], now: datetime | None = None) -> datetime:
    """Parse a `--delay-start` time with dateparser, preferring the future.

    Relative delays (`in 2 hours`, `90m`) and absolute times (`14:30`, `tomorrow 9am`,
    `2026-07-22 09:00`) are accepted; a time of day already past means tomorrow.

    Args:
        tokens: The option's words, joined so they need no quoting.
        now: Current time; defaults to now.

    Returns:
        The start time.

    Raises:
        ValueError: If the time cannot be parsed or is not in the future.
    """
    import dateparser

    text = " ".join(tokens).strip()
    if not text:
        raise ValueError("empty start time")
    now = now or datetime.now()
    when = dateparser.parse(
        text,
        settings={
            "PREFER_DATES_FROM": "future",
            "RELATIVE_BASE": now.replace(microsecond=0),
            "RETURN_AS_TIMEZONE_AWARE": False,
        },
    )
    if when is None:
        raise ValueError(
            f"could not understand start time {text!r} "
            f"(try e.g. 'in 2 hours', '90m', 'tomorrow 9am', '2026-07-22 09:00')"
        )
    if when <= now:
        raise ValueError(f"start time {when:%Y-%m-%d %H:%M:%S} is not in the future")
    return when


def _add_discover_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Add the `discover` subcommand's arguments."""
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to the config JSON (DiscoveryConfig schema).",
    )
    parser.add_argument(
        "--variant-density-blockfile",
        type=str,
        default="variant_density_blocks.tsv",
        help=(
            "Path to the variant density block file: "
            "whitespace-separated 3-column file CHROM START END (header optional)."
        ),
    )
    parser.add_argument(
        "--run-blockfile",
        type=str,
        default=None,
        help=(
            "Optional whitespace-separated 3-column file CHROM START END (header optional). "
            "If set, restricts step2 to within these blocks only (inclusive)."
        ),
    )
    parser.add_argument(
        "--start-phase",
        type=str,
        choices=PHASE_STRINGS_ORDERED + ["latest"],
        default="latest",
        help="Phase to start at (default: start or resume progress along --phases-to-run).",
    )
    parser.add_argument(
        "--phases-to-run",
        type=str,
        nargs="+",
        default=PHASE_STRINGS_ORDERED,
        choices=PHASE_STRINGS_ORDERED,
        help="Ordered space-separated list of phases to run (default: all phases).",
    )
    parser.add_argument(
        "--one-phase",
        action="store_true",
        help="Only run a single phase then exit.",
    )
    parser.add_argument(
        "--covar-stratum",
        type=str,
        nargs="+",
        default=None,
        help=(
            "With stratify_by_covar set in the config, restrict this invocation to the named "
            "covariate strata, 'combined' for the unsplit run (default: every covariate stratum "
            "the config asks for, run in parallel)."
        ),
    )
    parser.add_argument(
        "--meta-mode",
        type=str,
        choices=["auto"] + [m.value for m in MetaMode],
        default="auto",
        help=(
            "Specify the meta-analysis mode (default: auto). auto picks from --phases-to-run: "
            "singles_only if only singles phases (step2/conditioning) are listed, aggregates_only if "
            "only aggregates phases are, otherwise all."
        ),
    )
    parser.add_argument(
        "--overwrite-config",
        action="store_true",
        help="Overwrite persisted config values with those from the input --config.",
    )
    parser.add_argument(
        "--reset-jobs",
        type=int,
        default=None,
        help=(
            "State value to reset failed jobs to (default: previous state). "
            "In conditioning phases, acts as the conditioning JobState reset value."
        ),
    )
    parser.add_argument(
        "--force-reset",
        action="store_true",
        help=(
            "Force current phase into ERROR then apply --reset-jobs. "
            "In conditioning phases, forces the conditioning step into ERROR."
        ),
    )
    parser.add_argument(
        "--reset-jobs-list",
        nargs="+",
        type=str,
        default=None,
        help="Space-separated list of exact job descriptions to reset.",
    )
    parser.add_argument(
        "--force-job-retry",
        action="store_true",
        help=(
            "Resubmit failed jobs, with inputs rebuilt from the current config (in conditioning, a step "
            "that failed outside its jobs reruns the stage it failed in). Running and completed jobs "
            "are left alone, and if nothing failed the run just resumes."
        ),
    )
    parser.add_argument(
        "--list-jobs",
        action="store_true",
        help="List current jobs (descriptions and states) then exit.",
    )
    parser.add_argument(
        "--abort-all",
        action="store_true",
        help=(
            "Abort every running cloud job recorded in this run's state (all phases and conditioning "
            "steps), after confirmation, then exit. Stop any fedcat process driving the run first. "
            "State files are unchanged: on resume the aborted jobs are recorded as failed, so rerun "
            "them with --force-job-retry."
        ),
    )
    parser.add_argument(
        "--report",
        action="store_true",
        help=(
            "Report job states, task counts and costs for every phase and conditioning step of this run, "
            "then exit. Costs come from the platform where it reports them (DNAnexus) or from Cromwell's "
            "task metadata, which covers failed and preempted attempts too; where neither is reachable "
            "they are estimated from the logs of successful tasks."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Build jobs, serialise a state snapshot and print a summary, then exit without submitting.",
    )
    parser.add_argument(
        "--skip-validation",
        action="store_true",
        help=(
            "Skip the checks that workflows and remote input files (pheno, cov, pfile, step1 input, "
            "annotations) exist."
        ),
    )
    parser.add_argument(
        "--logging-level",
        type=str,
        default="INFO",
        choices=["NOFILE", "TRACE", "DEBUG", "INFO"],
        help="Logging verbosity (NOFILE → TRACE to stdout only, no log files).",
    )

    parser.add_argument(
        "--delay-start",
        nargs="+",
        default=None,
        metavar="WHEN",
        help=(
            "Defer job launch until a flexible time, given as a relative delay or an "
            "absolute time: 'in 2 hours', '90m', '1h30m', 'tomorrow 9am', '14:30', "
            "'2026-07-22 09:00' (see dateparser documentation)."
        ),
    )

    # ── Conditioning controls ─────────────────────────────────────────────────
    cond_grp = parser.add_argument_group(
        "Conditioning options",
        description=(
            "Runtime controls for the conditioning phases. "
            "Threshold parameters (p-value, r², collinearity checking) belong in "
            "the config JSON (cond_p_agg_retest, cond_p_singles_retest, "
            "cond_rsq_threshold, cond_check_coll, cond_max_steps)."
        ),
    )
    cond_grp.add_argument(
        "--cond-step",
        type=int,
        default=None,
        dest="cond_step",
        help="Start conditioning at this step number (default: 0 / autodetect from disk).",
    )
    cond_grp.add_argument(
        "--cond-reset-step",
        type=int,
        default=None,
        dest="cond_reset_step",
        help="Reset a failed conditioning step to this StepState integer value.",
    )
    cond_grp.add_argument(
        "--cond-one-step",
        action="store_true",
        dest="cond_one_step",
        help="Run only a single conditioning step then exit.",
    )
    cond_grp.add_argument(
        "--cond-singles-tarball",
        type=str,
        default=None,
        dest="cond_singles_tarball",
        metavar="FILE",
        help=(
            "Path to a pre-built singles subset results tarball to start step-0 singles "
            "conditioning from, instead of the one built by the step0_subset phase. "
            "Supported only for a single-phenotype run."
        ),
    )
    cond_grp.add_argument(
        "--cond-agg-tarball",
        type=str,
        default=None,
        dest="cond_agg_tarball",
        metavar="FILE",
        help=(
            "Path to a pre-built aggregate subset results tarball to start step-0 aggregate "
            "conditioning from, instead of the one built by the step0_subset phase. "
            "Supported only for a single-phenotype run."
        ),
    )
    return parser


def _log_run_summary(disc: DiscoveryRun, args: argparse.Namespace, requested_phase: Phase) -> None:
    """Log what the run will do (cohorts, phenotypes, phases, thresholds) and where its state is saved."""
    cfg = disc.config
    lines: list[str] = ["Run configuration summary:"]
    lines.append(f"    fedcat version : {__version__}")
    lines.append(f"    run_prefix     : {cfg.run_prefix}")
    if cfg.step1_run_prefix and cfg.step1_run_prefix != cfg.run_prefix:
        lines.append(f"    step1_run_prefix: {cfg.step1_run_prefix}  (reusing step1 output)")

    for bb, cohort in cfg.cohorts.items():
        strata = ", ".join(cohort.strata)
        sizes = ", ".join(f"{s}={cohort.cohort_sizes.get(s, '?')}" for s in cohort.strata)
        lines.append(f"    cohort {bb.value:<4}    : strata [{strata}]  sizes [{sizes}]")
        with contextlib.suppress(Exception):  # the summary must never break startup
            lines.append(f"                     output: {cfg.output_location(bb)}")

    phenos = cfg.pheno_names
    pheno_list = ", ".join(phenos) if len(phenos) <= 12 else ", ".join(phenos[:12]) + f", … (+{len(phenos) - 12})"
    lines.append(f"    phenotypes ({len(phenos)}){' ' if len(phenos) < 10 else ''}: {pheno_list}")
    lines.append(f"    chromosomes    : {', '.join(cfg.chromosomes_to_run)}")
    lines.append(
        f"    trait          : {cfg.trait_type.value} / {cfg.trait_transform.value}  regenie={cfg.regenie_version}  "
        f"coding_only={cfg.coding_only}  singles_mac>={cfg.singles_mac_threshold}  "
        f"aggregate_freqs={','.join(cfg.aggregate_freqs)}"
    )
    meta_location = cfg.meta_analysis_config.location.value if cfg.meta_analysis_config else "not configured"
    lines.append(f"    local env      : {cfg.local_env.value}")
    auto = " (auto)" if args.meta_mode == "auto" else ""
    lines.append(f"    meta loc/mode  : {meta_location} / {disc.meta_mode.value}{auto}")

    phase_str = " -> ".join(p.value for p in disc.phases_to_run)
    lines.append(
        f"    phases to run  : {phase_str}  (starting at: {requested_phase.value}"
        f"{', single phase only' if args.one_phase else ''})"
    )
    lines.append(
        f"    conditioning   : max_steps={cfg.cond_max_steps}  rsq<{cfg.cond_rsq_threshold}  "
        f"check_coll={cfg.cond_check_coll}  engine={cfg.cond_engine}"
        + (f" ({cfg.cond_mancojo_effect_mode})" if cfg.cond_engine == "mancojo" else "")
    )
    lines.append(
        f"                     p singles retest={cfg.cond_p_singles_retest:g} select={cfg.cond_p_singles_select:g}"
        f" | p agg retest={cfg.cond_p_agg_retest:g} select={cfg.cond_p_agg_select:g}"
    )
    if disc.cond_singles_tarball:
        lines.append(f"                     step0 singles tarball (direct): {disc.cond_singles_tarball}")
    if disc.cond_agg_tarball:
        lines.append(f"                     step0 agg tarball (direct)    : {disc.cond_agg_tarball}")

    if cfg.aou_bucket:
        lines.append(f"    aou_bucket     : gs://{cfg.aou_bucket.removeprefix('gs://')}")
    if BioBank.UKB in cfg.cohorts:
        lines.append(f"    dx_project_id  : {cfg.dx_project_id}")

    lines.append(f"    state dir      : {disc.run_path.resolve()}")
    cloud_state = disc._state_cloud_folder()
    lines.append(f"    cloud state    : {'gs://' + cloud_state if cloud_state else 'none (local only)'}")
    lines.append(
        f"    flags          : dry_run={args.dry_run} skip_validation={args.skip_validation} "
        f"overwrite_config={args.overwrite_config}"
    )

    logger.info("\n".join(lines))


def _run_context(config: DiscoveryConfig) -> str:
    """Return a run's logging context; covariate strata differ in `run_prefix`, so each has its own."""
    return f"{config.run_prefix}_discovery"


def _covar_stratum_label(config: DiscoveryConfig) -> str | None:
    """Return the covariate stratum shown in log line prefixes, or None for an unstratified run."""
    if config.stratify_by_covar is None:
        return None
    return config.covar_stratum or COVAR_STRATUM_COMBINED_NAME


def _select_covar_strata(
    config_obj: DiscoveryConfig, requested: list[str] | None
) -> list[tuple[str, DiscoveryConfig]]:
    """Return the `(covar_stratum, config)` pairs to run, restricted to `--covar-stratum` if given.

    The unstratified run is `combined` on the command line and `""` internally. Exits
    on a stratum the config does not define.
    """
    covar_strata = config_obj.covar_stratum_configs()
    if not requested:
        return covar_strata
    if config_obj.stratify_by_covar is None:
        logger.error("--covar-stratum requires stratify_by_covar to be set in the config.")
        sys.exit(1)
    available = [label or COVAR_STRATUM_COMBINED_NAME for label, _ in covar_strata]
    unknown = [s for s in requested if s not in available]
    if unknown:
        logger.error(
            f"Covariate strata {unknown} are not produced by this config; it runs: {', '.join(available)}"
        )
        sys.exit(1)
    wanted = {COVAR_STRATUM_COMBINED if s == COVAR_STRATUM_COMBINED_NAME else s for s in requested}
    return [(label, cfg) for label, cfg in covar_strata if label in wanted]


async def _discover_main(args: argparse.Namespace) -> None:
    """Run the `discover` subcommand: prepare each covariate stratum's run, then run them in parallel."""
    # ------------------------------------------------------------------
    # Load and validate the config
    # ------------------------------------------------------------------
    config_path = Path(args.config)
    if not config_path.exists():
        logger.error(f"Config file {config_path} not found.")
        sys.exit(1)

    conv = make_conv()
    raw_cfg = json.loads(config_path.read_text())
    try:
        config_obj: DiscoveryConfig = conv.structure(raw_cfg, DiscoveryConfig)
    except Exception:
        logger.exception(f"Failed to parse config file {config_path}")
        sys.exit(1)

    try:
        config_obj.validate()
    except Exception:
        logger.exception("Configuration validation failed")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Parse --delay-start now, to fail early; wait only after the checks pass
    # ------------------------------------------------------------------
    start_dt: datetime | None = None
    if args.delay_start is not None:
        try:
            start_dt = _parse_when(args.delay_start)
        except ValueError as exc:
            logger.error(f"Invalid --delay-start: {exc}")
            sys.exit(1)

    # ------------------------------------------------------------------
    # One run per covariate stratum, prepared in turn so prompts (phase
    # overwrite, bucket choice) come in order and at most once
    # ------------------------------------------------------------------
    covar_strata = _select_covar_strata(config_obj, args.covar_stratum)
    bucket_cache: dict[str, str] = {}
    prepared: list[DiscoveryRun] = []
    for covar_stratum, cfg in covar_strata:
        if len(covar_strata) > 1:
            logger.info(f"── Preparing {covar_stratum or COVAR_STRATUM_COMBINED_NAME} run ({cfg.run_prefix}) ──")
        disc = await _prepare_discovery_run(args, cfg, conv, bucket_cache)
        if disc is not None:
            prepared.append(disc)

    if not prepared:
        return

    # ------------------------------------------------------------------
    # Wait for --delay-start
    # ------------------------------------------------------------------
    if start_dt is not None:
        remaining = start_dt - datetime.now()
        if remaining.total_seconds() <= 0:
            logger.info(f"Scheduled start {start_dt:%Y-%m-%d %H:%M:%S} is not in the future; launching now.")
        else:
            logger.info(
                f"Waiting until {start_dt:%Y-%m-%d %H:%M:%S} to launch jobs "
                f"({timedelta(seconds=int(remaining.total_seconds()))} from now)."
            )
            await asyncio.sleep(remaining.total_seconds())
            logger.info("Scheduled start time reached; launching jobs.")

    # ------------------------------------------------------------------
    # Run concurrently; one failing does not cancel the others
    # ------------------------------------------------------------------
    async def _run_one(disc: DiscoveryRun) -> None:
        """Run one covariate stratum's run, logging in its own context."""
        # Set in the task, whose context copy everything it starts inherits
        setup_run_logging(run_context=_run_context(disc.config), covar_stratum=_covar_stratum_label(disc.config))
        await disc.run(
            one_phase=args.one_phase,
            reset_jobs=JobState(args.reset_jobs) if args.reset_jobs is not None else None,
            reset_jobs_list=args.reset_jobs_list,
            force_reset=args.force_reset,
        )

    failed = False
    STATE_WRITE_QUEUE.start()
    try:
        if len(prepared) == 1:
            await _run_one(prepared[0])
        else:
            results = await asyncio.gather(*(_run_one(d) for d in prepared), return_exceptions=True)
            for disc, result in zip(prepared, results):
                if isinstance(result, BaseException):
                    logger.error(f"Run {disc.config.run_prefix!r} raised: {result}")
                    failed = True
    finally:
        await STATE_WRITE_QUEUE.stop()

    await logger.complete()
    # A non-zero exit lets scripts that chain runs stop on a failed phase
    if failed or any(disc.had_error for disc in prepared):
        sys.exit(1)


async def _prepare_discovery_run(
    args: argparse.Namespace,
    config_obj: DiscoveryConfig,
    conv,
    bucket_cache: dict[str, str],
) -> DiscoveryRun | None:
    """Prepare one discovery run: phases, workspace bucket, checks, references and state.

    `--list-jobs`, `--dry-run`, `--abort-all` and `--report` are handled here. Exits on
    errors.

    Args:
        args: Parsed `discover` arguments.
        config_obj: The covariate stratum's config.
        conv: Converter from [`make_conv`][fedcat.core.structures.make_conv].
        bucket_cache: Workspace bucket chosen for an earlier covariate stratum of this
            invocation, shared so it is asked for at most once.

    Returns:
        The run, ready to run, or None when there is nothing to run.
    """
    # ------------------------------------------------------------------
    # Run directory and logging
    # ------------------------------------------------------------------
    run_dir = config_obj.discovery_run_path
    run_dir.mkdir(parents=True, exist_ok=True)

    if args.logging_level == "NOFILE":
        effective_level = "TRACE"
        base_log_path = None
    else:
        effective_level = args.logging_level
        base_log_path = run_dir / "logs"

    setup_run_logging(
        run_context=_run_context(config_obj),
        base_log_path=base_log_path,
        log_warnings_separately=False,
        logging_level=effective_level,
        covar_stratum=_covar_stratum_label(config_obj),
    )
    STAGE_VAR.set("setup")

    logger.info(f"Starting fedcat version {__version__}")
    logger.debug(f"Run directory: {run_dir.resolve()}")

    # ------------------------------------------------------------------
    # Phases to run, and the phase to start from
    # ------------------------------------------------------------------
    phases_to_run = [Phase(p) for p in args.phases_to_run]

    phases_to_run_ordered = [p for p in PHASE_ORDERED if p in phases_to_run]
    if phases_to_run != phases_to_run_ordered:
        logger.error(
            f"--phases-to-run {args.phases_to_run} not in valid order; "
            f"try: {' '.join(p.value for p in phases_to_run_ordered)}"
        )
        sys.exit(1)

    meta_mode = MetaMode.from_phases(phases_to_run) if args.meta_mode == "auto" else MetaMode(args.meta_mode)

    # An interaction run reports several TEST rows per variant (ADD-CONDTL, ADD-INT_*).
    # The meta-analysis handles them, but the step-0 subset and conditioning do not
    if config_obj.interaction:
        downstream = [p for p in phases_to_run if p in (Phase.COND_SINGLES, Phase.COND_AGG)]
        if downstream:
            logger.error(
                f"interaction={config_obj.interaction!r} is set and "
                f"{', '.join(p.value for p in downstream)} requested: "
                "conditioning is not currently supported for interaction runs."
            )
            sys.exit(1)
        if Phase.STEP2_AGG in phases_to_run:
            logger.error(
                f"interaction={config_obj.interaction!r} is set and {Phase.STEP2_AGG.value} requested: "
                "aggregates are not currently supported for interaction runs."
            )
            sys.exit(1)

    existing_phases: list[Phase] = [
        p for p in phases_to_run if (run_dir / f"discovery_state_{p.value}.json").exists()
    ]

    if args.start_phase == "latest":
        if existing_phases:
            requested_phase = existing_phases[-1]
            logger.info(f"Resuming from latest existing phase: {requested_phase.value}")
        else:
            requested_phase = phases_to_run[0]
            logger.info(
                f"No existing state files found; starting at first requested phase: {requested_phase.value}"
            )
    else:
        requested_phase = Phase(args.start_phase)
        if requested_phase not in phases_to_run:
            logger.error(
                f"--start-phase {requested_phase.value} not in --phases-to-run {[p.value for p in phases_to_run]}"
            )
            sys.exit(1)
        elif existing_phases and (
            PHASE_ORDERED.index(requested_phase) < PHASE_ORDERED.index(existing_phases[-1])
        ):
            logger.warning(
                f"--start-phase {requested_phase.value} is before the latest existing phase "
                f"{existing_phases[-1].value}."
            )
            answer = input(f"Continue from {requested_phase.value} and overwrite later phases? y/N: ")
            if answer.strip().lower() not in ("y", "yes"):
                logger.info("Exiting.")
                sys.exit(0)
            _set_aside_later_phases(run_dir, requested_phase)
            logger.info(f"Continuing from {requested_phase.value}; later phases will be overwritten.")
        else:
            logger.info(f"Starting at requested phase: {requested_phase.value}")

    if requested_phase in (Phase.COND_SINGLES, Phase.COND_AGG) and args.reset_jobs_list is not None:
        logger.error("--reset-jobs-list not implemented for conditioning phases")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Workspace bucket: from the config, the saved state, or chosen by the user.
    # A configured bucket is looked up with `wb` for its resource id.
    # ------------------------------------------------------------------
    if not args.list_jobs and BioBank.AOU in config_obj.cohorts and not config_obj.aou_bucket:
        if cached := bucket_cache.get("bucket"):
            config_obj.aou_bucket = cached
            config_obj.aou_bucket_resource_id = bucket_cache.get("resource_id", "")
            logger.info(f"Reusing workspace bucket resolved for this invocation: {cached}")
        elif persisted := _persisted_workspace_bucket(run_dir, existing_phases):
            config_obj.aou_bucket, config_obj.aou_bucket_resource_id = persisted
            logger.info(
                f"Reusing saved workspace bucket: {config_obj.aou_bucket} "
                f"(resource id: {config_obj.aou_bucket_resource_id})"
            )
        else:
            from fedcat.platforms.gcs_cromwell.buckets import resolve_workspace_bucket

            logger.info("Searching for available workspace buckets...")
            chosen = resolve_workspace_bucket()
            config_obj.aou_bucket = chosen.gs_url
            config_obj.aou_bucket_resource_id = chosen.resource_id
            logger.info(f"Using workspace bucket: {chosen.gs_url} (resource id: {chosen.resource_id})")
        bucket_cache.setdefault("bucket", config_obj.aou_bucket)
        bucket_cache.setdefault("resource_id", config_obj.aou_bucket_resource_id)
    elif (
        not args.list_jobs
        and BioBank.AOU in config_obj.cohorts
        and config_obj.cohorts[BioBank.AOU].workflow_backend is WorkflowBackend.WB_CLI
    ):
        from fedcat.platforms.gcs_cromwell.buckets import resolve_configured_workspace_bucket

        configured = resolve_configured_workspace_bucket(
            config_obj.aou_bucket,
            config_obj.aou_bucket_resource_id,
        )
        config_obj.aou_bucket_resource_id = configured.resource_id
        logger.info(
            f"Using configured workspace bucket: {config_obj.aou_bucket} "
            f"(resource id: {config_obj.aou_bucket_resource_id})"
        )

    # ------------------------------------------------------------------
    # --abort-all and --report need only the state files and the platforms
    # ------------------------------------------------------------------
    if args.abort_all:
        await _abort_all(config_obj, run_dir, conv)
        return None
    if args.report:
        await _report(config_obj, run_dir, conv)
        return None

    # ------------------------------------------------------------------
    # Preflight checks, also for a dry run. They connect the platforms, which
    # resolving references below needs.
    # ------------------------------------------------------------------
    if not args.list_jobs:
        from fedcat.core.checks import (
            CheckContext,
            RunFlags,
            log_findings,
            run_checks,
        )
        from fedcat.platforms import build_platforms

        blockfile = Path(args.variant_density_blockfile) if args.variant_density_blockfile else None
        ctx = CheckContext(
            config=config_obj,
            platforms=build_platforms(config_obj),
            phases_to_run=phases_to_run,
            variant_density_blockfile=blockfile,
            meta_mode=meta_mode,
            flags=RunFlags(dry_run=args.dry_run, skip_validation=args.skip_validation),
        )
        if log_findings(await run_checks(ctx)):
            sys.exit(1)

    # ------------------------------------------------------------------
    # Put each cohort's references in the form its platform's jobs take:
    # files as data refs, templates and annotations as path refs
    # ------------------------------------------------------------------
    if not args.list_jobs:
        from fedcat.platforms import build_platforms

        platforms = build_platforms(config_obj)
        for bb, cohort in config_obj.cohorts.items():
            platform = platforms[bb]
            cohort.pheno_file = platform.canonical_data_ref(cohort.pheno_file)
            cohort.cov_file = platform.canonical_data_ref(cohort.cov_file)
            if cohort.cond_file:
                cohort.cond_file = platform.canonical_data_ref(cohort.cond_file)
            cohort.pfile_template = platform.canonical_path_ref(cohort.pfile_template)
            cohort.step1_input_template = platform.canonical_path_ref(cohort.step1_input_template)
            for annots in cohort.regenie_annotations:
                annots.mask_file = platform.canonical_path_ref(annots.mask_file)
                annots.annot_file = platform.canonical_path_ref(annots.annot_file)
                annots.setlist_file = platform.canonical_path_ref(annots.setlist_file)

    if not args.list_jobs and (metaconf := config_obj.meta_analysis_config) is not None:
        match metaconf.location:
            case RunEnv.AOU:
                if metaconf.schema_directory:
                    metaconf.schema_directory = "gs://" + metaconf.schema_directory.removeprefix("gs://")
                if metaconf.regenie_annotations:
                    for annots in metaconf.regenie_annotations:
                        annots.mask_file = "gs://" + annots.mask_file.removeprefix("gs://")
                        annots.annot_file = "gs://" + annots.annot_file.removeprefix("gs://")
                        annots.setlist_file = "gs://" + annots.setlist_file.removeprefix("gs://")
            case RunEnv.LOCAL:
                pass

    # ------------------------------------------------------------------
    # Load the phase's state, or create the run
    # ------------------------------------------------------------------
    state_path = run_dir / f"discovery_state_{requested_phase.value}.json"

    if args.list_jobs:
        if not state_path.exists():
            logger.error(f"No existing state file {state_path} to list jobs from.")
            return None
        disc: DiscoveryRun = _load_discovery_state(state_path, conv)
        logger.info(f"Total of {len(disc.jobs)} jobs for phase {disc.current_phase.value}")
        for i, job in enumerate(disc.jobs):
            logger.info(f"JOB {i + 1}: {job.description}, state={job.state.name}")
        await logger.complete()
        return None

    if state_path.exists():
        disc: DiscoveryRun = _load_discovery_state(state_path, conv)
        if disc.current_phase != requested_phase:
            logger.error(
                f"{state_path} holds {disc.current_phase.value} state, not {requested_phase.value}; the file "
                f"does not match its name, so fix or remove it rather than resume from it."
            )
            sys.exit(1)
        from fedcat.core.engine import (
            adopt_run_controls,
            config_differences,
            format_config_differences,
        )

        diffs = config_differences(config_obj, disc.config)
        adopt_run_controls(disc.config, config_obj)
        if args.overwrite_config:
            disc.adopt_config(config_obj)
            if diffs:
                logger.info(f"Updated the saved config from {args.config}:\n{format_config_differences(diffs)}")
            else:
                logger.info(f"--overwrite-config: {args.config} matches the saved config, nothing to update")
        elif diffs:
            logger.error(
                f"The run was saved with a different config from {args.config}:\n"
                f"{format_config_differences(diffs)}\nUse --overwrite-config to replace it."
            )
            sys.exit(1)
        disc.phases_to_run = phases_to_run
        logger.info(f"Resumed discovery state from {state_path}")
    else:
        disc = DiscoveryRun(
            config=config_obj,
            run_path=run_dir,
            variant_density_blockfile=Path(args.variant_density_blockfile),
            phases_to_run=phases_to_run,
            current_phase=requested_phase,
        )

    disc.cond_step = getattr(args, "cond_step", None)
    disc.cond_one_step = getattr(args, "cond_one_step", False)
    disc.cond_reset_step = getattr(args, "cond_reset_step", None)
    disc._overwrite_config = args.overwrite_config

    for flag, dest in (
        ("--cond-singles-tarball", "cond_singles_tarball"),
        ("--cond-agg-tarball", "cond_agg_tarball"),
    ):
        if (value := getattr(args, dest, None)) is not None:
            tarball = Path(value)
            if not tarball.exists():
                logger.error(f"{flag} file {tarball} not found.")
                sys.exit(1)
            setattr(disc, dest, tarball.resolve())

    # In conditioning phases the reset flags apply to the conditioning steps
    if requested_phase in (Phase.COND_SINGLES, Phase.COND_AGG):
        if args.force_job_retry:
            disc.cond_retry_failed = True
        else:
            disc.cond_force_reset = args.force_reset
            disc.cond_reset_jobs = args.reset_jobs

    # ------------------------------------------------------------------
    # Options applied whether resuming or not
    # ------------------------------------------------------------------
    if args.run_blockfile:
        if not Path(args.run_blockfile).exists():
            logger.error(f"Blockfile {args.run_blockfile} not found.")
            sys.exit(1)
        else:
            disc.run_blockfile = Path(args.run_blockfile)
            logger.info(f"Using blockfile {args.run_blockfile} to restrict step2 regions.")

    disc.meta_mode = meta_mode

    _log_run_summary(disc, args, requested_phase)

    # ------------------------------------------------------------------
    # Dry run: build the jobs and save a snapshot
    # ------------------------------------------------------------------
    if args.dry_run:
        disc.jobs = await disc._build_jobs_for_phase()
        await disc._snapshot_state("dry_run_no_overwrite")
        preview = 10
        logger.info(f"Total of {len(disc.jobs)} jobs for phase {disc.current_phase.value}")
        for i, job in enumerate(disc.jobs[:preview]):
            logger.info(f"JOB {i + 1}: {job.description}, state={job.state.name}")
        if len(disc.jobs) > preview:
            logger.info(f"... ({len(disc.jobs) - preview} more jobs not shown)")
        await logger.complete()
        return None

    # Make the local input files readable
    local_prereqs: list[str] = []
    if Phase.STEP2_SINGLES in disc.phases_to_run or Phase.STEP2_AGG in disc.phases_to_run:
        local_prereqs.append(args.variant_density_blockfile)
    for file in local_prereqs:
        p = Path(file)
        if p.exists():
            p.chmod(0o744)

    return disc


def _convert_annots_main(args: argparse.Namespace) -> None:
    """Run the `convert-annots` subcommand."""
    from fedcat_rs import convert_tsv_to_parquet_rs

    masks_file = [Path(m) for m in args.masks]
    annots_files = [Path(a) for a in args.annots]
    output_dir = Path(args.output_dir)

    missing = [str(f) for f in annots_files if not f.exists()] + [str(m) for m in masks_file if not m.exists()]
    if missing:
        logger.error(f"File(s) not found: {missing}")
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    annotation_tsvs = [str(f) for f in annots_files]
    mask_tsvs = [str(m) for m in masks_file]
    logger.info(f"Converting {len(annotation_tsvs)} annotation TSV file(s) to parquet schema in {output_dir}")

    meta = convert_tsv_to_parquet_rs(annotation_tsvs, mask_tsvs, str(output_dir), None)
    logger.info(
        f"Conversion complete: {meta['n_rows']} rows, "
        f"{meta['n_variants']} variants, "
        f"{meta['n_units']} units, "
        f"{meta['n_annots']} annotations"
    )
    logger.info(f"Schema written to: {output_dir.resolve()}")


def _add_convert_annots_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Add the `convert-annots` subcommand's arguments."""
    parser.add_argument(
        "--masks",
        type=str,
        nargs="+",
        required=True,
        metavar="FILE",
        help="One or more REGENIE masks definition TSV files (MASK<tab>annot1,annot2,...).",
    )
    parser.add_argument(
        "--annots",
        type=str,
        nargs="+",
        required=True,
        metavar="FILE",
        help="One or more REGENIE annotation TSV files (VARIANT<tab>UNIT<tab>ANNOT).",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        metavar="DIR",
        dest="output_dir",
        help="Directory to write the parquet schema files into (created if absent).",
    )
    return parser


def _add_meta_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Add the `meta` subcommand's arguments."""
    parser.add_argument(
        "--cohort-file",
        type=str,
        required=True,
        metavar="FILE",
        dest="cohort_file",
        help=(
            "Text file defining per-cohort input files. "
            "Lines starting with '# <name>' define cohort names; "
            "subsequent non-empty lines are file paths for that cohort."
        ),
    )
    parser.add_argument(
        "--output-file",
        type=str,
        required=True,
        metavar="FILE",
        dest="output_file",
        help="Output file path for meta-analysis results.",
    )
    parser.add_argument(
        "--res-type",
        type=str,
        required=True,
        metavar="TYPE",
        dest="res_type",
        help="Result type (e.g. 'single_variants', 'burden', 'mac10').",
    )
    parser.add_argument(
        "--compression",
        type=str,
        default=None,
        choices=["zst", "gz"],
        help="Output compression codec. If omitted, inferred from output file extension.",
    )
    parser.add_argument(
        "--compression-level",
        type=int,
        default=None,
        metavar="N",
        dest="compression_level",
        help="Compression level (codec-specific).",
    )
    parser.add_argument(
        "--backend",
        type=str,
        default="streaming",
        choices=["streaming", "merge"],
        help="Meta-analysis backend (default: streaming).",
    )
    parser.add_argument(
        "--aaf-bins",
        type=str,
        default=None,
        metavar="F[,F...]",
        dest="aaf_bins",
        help=(
            "Comma-separated allele-frequency bins. Each cohort contributes to bin F only "
            "where its own A1FREQ <= F; output has one row per variant and bin (AAF_BIN "
            "column). For the single-variant input to ACAT."
        ),
    )
    parser.add_argument("--chrom-col", type=str, default=None, dest="col_chrom", help="Column name for CHROM.")
    parser.add_argument("--genpos-col", type=str, default=None, dest="col_genpos", help="Column name for GENPOS.")
    parser.add_argument("--id-col", type=str, default=None, dest="col_id", help="Column name for ID.")
    parser.add_argument(
        "--allele0-col", type=str, default=None, dest="col_allele0", help="Column name for ALLELE0."
    )
    parser.add_argument(
        "--allele1-col", type=str, default=None, dest="col_allele1", help="Column name for ALLELE1."
    )
    parser.add_argument("--a1freq-col", type=str, default=None, dest="col_a1freq", help="Column name for A1FREQ.")
    parser.add_argument("--n-col", type=str, default=None, dest="col_n", help="Column name for N.")
    parser.add_argument("--beta-col", type=str, default=None, dest="col_beta", help="Column name for BETA.")
    parser.add_argument("--se-col", type=str, default=None, dest="col_se", help="Column name for SE.")
    parser.add_argument("--log10p-col", type=str, default=None, dest="col_log10p", help="Column name for LOG10P.")
    parser.add_argument("--test-col", type=str, default=None, dest="col_test", help="Column name for TEST.")
    parser.add_argument(
        "--col-separator",
        type=str,
        default=None,
        dest="col_separator",
        help="Field separator character for input files.",
    )
    return parser


def _meta_main(args: argparse.Namespace) -> None:
    """Run the `meta` subcommand."""
    from fedcat.core.io import parse_cohort_file
    from fedcat.stats.meta_analysis import meta_analyse_regenie

    cohort_path = Path(args.cohort_file)
    if not cohort_path.exists():
        logger.error(f"Cohort file not found: {cohort_path}")
        sys.exit(1)

    per_cohort_files = parse_cohort_file(cohort_path)
    output_file = Path(args.output_file)
    # Without inputs (e.g. no mac10 results) no output is written, which the workflows check for
    if not any(per_cohort_files.values()):
        logger.warning(f"No {args.res_type} files listed in {cohort_path}; writing no output")
        return

    col_keys = [
        ("chrom", args.col_chrom),
        ("genpos", args.col_genpos),
        ("id", args.col_id),
        ("allele0", args.col_allele0),
        ("allele1", args.col_allele1),
        ("a1freq", args.col_a1freq),
        ("n", args.col_n),
        ("beta", args.col_beta),
        ("se", args.col_se),
        ("log10p", args.col_log10p),
        ("test", args.col_test),
        ("separator", args.col_separator),
    ]
    column_config = {k: v for k, v in col_keys if v is not None} or None

    meta_analyse_regenie(
        per_cohort_files=per_cohort_files,
        output_file=output_file,
        res_type=args.res_type,
        column_config=column_config,
        compression=args.compression,
        compression_level=args.compression_level,
        backend=args.backend,
        aaf_bins=[f.strip() for f in args.aaf_bins.split(",") if f.strip()] if args.aaf_bins else None,
    )


def _add_acat_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Add the `acat` subcommand's arguments."""
    parser.add_argument(
        "--single-variants-file",
        type=str,
        required=True,
        metavar="FILE",
        dest="single_variants_file",
        help="TSV file with single-variant summary stats.",
    )
    parser.add_argument(
        "--burden-file",
        type=str,
        default=None,
        metavar="FILE",
        dest="burden_file",
        help="Optional TSV file with burden summary stats.",
    )
    parser.add_argument(
        "--output-file",
        type=str,
        required=True,
        metavar="FILE",
        dest="output_file",
        help="Output file path for ACAT results.",
    )
    # Annotations: --schema-dir, or --masks-tsvs with --annotations-tsvs
    parser.add_argument(
        "--schema-dir",
        type=str,
        default=None,
        metavar="DIR",
        dest="schema_dir",
        help="Directory containing pre-converted parquet schema (mutually exclusive with --masks-tsvs/--annotations-tsvs).",
    )
    parser.add_argument(
        "--masks-tsvs",
        type=str,
        nargs="+",
        default=None,
        metavar="FILE",
        dest="masks_tsvs",
        help="One or more REGENIE masks definition TSV files.",
    )
    parser.add_argument(
        "--annotations-tsvs",
        type=str,
        nargs="+",
        default=None,
        metavar="FILE",
        dest="annotations_tsvs",
        help="One or more REGENIE annotation TSV files.",
    )
    parser.add_argument(
        "--aggregate-freqs",
        type=str,
        default="0.001",
        metavar="F[,F...]",
        dest="aggregate_freqs",
        help=(
            "Comma-separated upper allele-frequency bounds for aggregate tests, used "
            "verbatim as threshold/ID labels (default: '0.001'). One set of ACAT "
            "aggregates is produced per threshold."
        ),
    )
    parser.add_argument("--id-col", type=str, default="ID", dest="id_col", help="ID column name (default: ID).")
    parser.add_argument(
        "--allele-freq-col",
        type=str,
        default="A1FREQ",
        dest="allele_freq_col",
        help="Allele frequency column name (default: A1FREQ).",
    )
    parser.add_argument(
        "--log10p-col",
        type=str,
        default="LOG10P",
        dest="log10p_col",
        help="-log10(P) column name (default: LOG10P).",
    )
    parser.add_argument(
        "--n-threads",
        type=int,
        default=None,
        metavar="N",
        dest="n_threads",
        help="Number of threads (default: cpu_count - 1).",
    )
    parser.add_argument(
        "--weighting-scheme",
        type=str,
        default="standard",
        choices=["standard", "equal"],
        dest="weighting_scheme",
        help="Variant weighting scheme (default: standard).",
    )
    parser.add_argument(
        "--compression",
        type=str,
        default=None,
        choices=["zst", "gz"],
        help="Output compression codec. If omitted, inferred from output file extension.",
    )
    parser.add_argument(
        "--compression-level",
        type=int,
        default=None,
        metavar="N",
        dest="compression_level",
        help="Compression level (codec-specific).",
    )
    parser.add_argument(
        "--backend",
        type=str,
        default="sweep",
        choices=["sweep", "arrow"],
        help="ACAT backend (default: sweep).",
    )
    parser.add_argument(
        "--no-full-bounds",
        action="store_false",
        dest="full_bounds",
        help=(
            "With TSV annotations, take mask START/END only from variants with single-variant results "
            "(faster; burden-only masks get NA). Default: bounds from all annotated variants."
        ),
    )
    return parser


def _acat_main(args: argparse.Namespace) -> None:
    """Run the `acat` subcommand."""
    from fedcat.stats.acat import compute_acat

    compute_acat(
        single_variants_file=Path(args.single_variants_file),
        burden_file=Path(args.burden_file) if args.burden_file else None,
        schema_dir=Path(args.schema_dir) if args.schema_dir else None,
        masks_tsvs=[Path(m) for m in args.masks_tsvs] if args.masks_tsvs else None,
        annotations_tsvs=[Path(a) for a in args.annotations_tsvs] if args.annotations_tsvs else None,
        output_file=Path(args.output_file),
        aggregate_freqs=[f.strip() for f in args.aggregate_freqs.split(",") if f.strip()],
        id_col=args.id_col,
        allele_freq_col=args.allele_freq_col,
        log10p_col=args.log10p_col,
        n_threads=args.n_threads,
        weighting_scheme=args.weighting_scheme,
        compression=args.compression,
        compression_level=args.compression_level,
        backend=args.backend,
        full_bounds=args.full_bounds,
    )


def _add_subset_results_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Add the `subset-results` subcommand's arguments."""
    parser.add_argument(
        "--phenotypes",
        type=str,
        nargs="+",
        required=True,
        metavar="PHENO",
        help="Phenotype names to process.",
    )
    parser.add_argument(
        "--chr",
        type=str,
        default=None,
        dest="chr_",
        metavar="CHR",
        help="Optional chromosome label for tarball naming (e.g. 1, 22, X).",
    )
    # The result files given decide which subsets are written
    parser.add_argument(
        "--acat-file",
        type=str,
        default=None,
        dest="acat_file",
        metavar="FILE",
        help="ACAT output TSV. If provided, ACAT hits are included in the aggregate tarball.",
    )
    parser.add_argument(
        "--burden-file",
        type=str,
        default=None,
        dest="burden_file",
        metavar="FILE",
        help="Burden meta-analysis TSV. If provided, burden hits are included in the aggregate tarball.",
    )
    parser.add_argument(
        "--singles-file",
        type=str,
        default=None,
        dest="singles_file",
        metavar="FILE",
        help="Unconditional singles meta-analysis TSV. If provided, singles subset tarball is created.",
    )
    # Annotations: --schema-dir, or --annotations with --masks-tsvs
    parser.add_argument(
        "--annotations",
        type=str,
        nargs="+",
        default=[],
        metavar="FILE",
        help="Annotation TSV files (required for --acat-file/--burden-file unless --schema-dir is given).",
    )
    parser.add_argument(
        "--masks-tsvs",
        type=str,
        nargs="+",
        default=[],
        dest="masks_tsvs",
        metavar="FILE",
        help="REGENIE masks TSV files (used with --annotations; ignored when --schema-dir is set).",
    )
    parser.add_argument(
        "--schema-dir",
        type=str,
        default=None,
        dest="schema_dir",
        metavar="DIR",
        help="Parquet schema directory (alternative to --annotations).",
    )
    parser.add_argument(
        "--log10p-threshold",
        type=float,
        default=6.0,
        metavar="T",
        dest="log10p_threshold",
        help="Minimum -log10(P) to include in tophits (default: 6.0 = P<1e-6).",
    )
    parser.add_argument(
        "--agg-p-threshold",
        type=float,
        default=None,
        metavar="P",
        dest="agg_p_threshold",
        help="P-value threshold for ACAT/burden tophits, i.e. cond_p_agg_retest (overrides --log10p-threshold).",
    )
    parser.add_argument(
        "--singles-p-threshold",
        type=float,
        default=None,
        metavar="P",
        dest="singles_p_threshold",
        help="P-value threshold for singles tophits, i.e. cond_p_singles_retest (overrides --log10p-threshold).",
    )
    return parser


def _subset_results_main(args: argparse.Namespace) -> None:
    """Run the `subset-results` subcommand."""
    from fedcat.stats.subset_results import subset_results_for_phenotype

    annotation_paths = [Path(a) for a in args.annotations] if args.annotations else None
    masks_paths = [Path(m) for m in args.masks_tsvs] if args.masks_tsvs else None
    schema_dir = Path(args.schema_dir) if args.schema_dir else None
    acat_file = Path(args.acat_file) if args.acat_file else None
    burden_file = Path(args.burden_file) if args.burden_file else None
    singles_file = Path(args.singles_file) if args.singles_file else None

    # Only the aggregate subset needs annotations
    if (acat_file or burden_file) and not annotation_paths and not schema_dir:
        logger.error("At least one of --annotations or --schema-dir is required with --acat-file/--burden-file")
        sys.exit(1)

    def _log10p(p: float | None) -> float:
        """Return -log10 of a p-value threshold, or `--log10p-threshold` if None; exits if not in (0, 1]."""
        if p is None:
            return args.log10p_threshold
        if not 0 < p <= 1:
            logger.error(f"P-value thresholds must be in (0, 1], got {p}")
            sys.exit(1)
        return -math.log10(p)

    agg_threshold = _log10p(args.agg_p_threshold)
    singles_threshold = _log10p(args.singles_p_threshold)

    for pheno in args.phenotypes:
        subset_results_for_phenotype(
            pheno,
            acat_path=acat_file,
            burden_path=burden_file,
            singles_path=singles_file,
            annotation_paths=annotation_paths,
            masks_paths=masks_paths,
            schema_dir=schema_dir,
            agg_threshold=agg_threshold,
            singles_threshold=singles_threshold,
            chr_label=args.chr_,
        )

    chr_msg = f" for chr{args.chr_}" if args.chr_ else ""
    logger.success(f"Subsetting complete{chr_msg} ({len(args.phenotypes)} phenotype(s))")


def _add_summarise_indep_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Add the `summarise-indep` subcommand's arguments."""
    parser.add_argument(
        "--run-dir",
        type=str,
        required=True,
        metavar="DIR",
        dest="run_dir",
        help="Conditioning run directory containing the {cond}_cond_step{N} subdirectories.",
    )
    parser.add_argument(
        "--cond-type",
        type=str,
        required=True,
        choices=["singles", "agg"],
        dest="cond_type",
        help="Conditioning type: 'singles' (single-variant) or 'agg' (aggregate).",
    )
    parser.add_argument(
        "--pheno",
        type=str,
        required=True,
        metavar="NAME",
        help="Phenotype name (used in the output filenames).",
    )
    parser.add_argument(
        "--passing-ids-json",
        type=str,
        default=None,
        metavar="FILE",
        dest="passing_ids_json",
        help="Override path to the JSON containing top-level 'passing_ids'. "
        "Auto-located within --run-dir when omitted.",
    )
    parser.add_argument(
        "--end-step",
        type=int,
        default=None,
        metavar="N",
        dest="end_step",
        help="Override the highest conditioning step to include. "
        "Auto-detected from the step directories when omitted.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        metavar="DIR",
        dest="output_dir",
        help="Override the results output directory (default: {run-dir}/{cond}_cond_results).",
    )
    parser.add_argument(
        "--stat-col",
        type=str,
        default=None,
        choices=["LOG10P", "P"],
        dest="stat_col",
        help="Test-statistic column base name. Auto-detected from the summary_output "
        "header when omitted (older runs use 'P' instead of 'LOG10P').",
    )
    return parser


def _summarise_indep_main(args: argparse.Namespace) -> None:
    """Run the `summarise-indep` subcommand."""
    from fedcat.core.structures import CondType
    from fedcat.stats.summarise_indep import summarise_independent_results

    run_dir = Path(args.run_dir)
    if not run_dir.is_dir():
        logger.error(f"Run directory not found: {run_dir}")
        sys.exit(1)

    passing_ids_json = Path(args.passing_ids_json) if args.passing_ids_json else None
    output_dir = Path(args.output_dir) if args.output_dir else None

    try:
        summarise_independent_results(
            run_dir=run_dir,
            cond_type=CondType(args.cond_type),
            pheno_name=args.pheno,
            passing_ids_json=passing_ids_json,
            end_step=args.end_step,
            output_dir=output_dir,
            stat_col=args.stat_col,
        )
    except (FileNotFoundError, ValueError) as e:
        logger.error(str(e))
        sys.exit(1)


def _add_collinearity_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Add the `collinearity` subcommand's arguments."""
    parser.add_argument("-p", "--pfile", type=str, required=True, help="plink2 pgen/pvar prefix.")
    parser.add_argument(
        "-s", "--selected", type=str, required=True, help="File of selected (conditioned-on) variant IDs."
    )
    parser.add_argument("-c", "--candidate", type=str, required=True, help="File of candidate variant IDs.")
    parser.add_argument("-o", "--out", type=str, required=True, help="Output TSV path.")
    parser.add_argument("-t", "--threads", type=int, required=True, help="BLAS thread count.")
    parser.add_argument(
        "-b",
        "--batch-size",
        type=int,
        default=1000,
        dest="batch_size",
        help="Candidate variants processed per batch (default: 1000).",
    )
    parser.add_argument(
        "--rsq-threshold",
        type=float,
        default=0.99,
        dest="rsq_threshold",
        help="Pairwise r² above which a candidate is flagged degenerate (default: 0.99).",
    )
    return parser


def _collinearity_main(args: argparse.Namespace) -> None:
    """Run the `collinearity` subcommand; BLAS threads are set from `--threads` on import."""
    from fedcat.stats.collinearity import run_collinearity_check

    run_collinearity_check(
        pfile=args.pfile,
        selected=args.selected,
        candidate=args.candidate,
        out=args.out,
        batch_size=args.batch_size,
        rsq_threshold=args.rsq_threshold,
    )


def _add_setup_workflows_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Add the `setup-workflows` subcommand's arguments."""
    parser.add_argument(
        "--biobank",
        type=str,
        default=BioBank.AOU.value,
        choices=[b.value for b in BioBank],
        help="Biobank workflow set to stage (default: aou). Only 'aou' is implemented.",
    )
    parser.add_argument(
        "--bucket",
        type=str,
        default=None,
        help="Workspace bucket (gs://… or bare name) for uploads. "
        "If omitted, you'll choose from the workspace's CONTROLLED buckets (via wb).",
    )
    parser.add_argument(
        "--bucket-id",
        type=str,
        default=None,
        dest="bucket_id",
        help="Workbench bucket resource id for 'wb workflow create'. "
        "Defaults to the resolved bucket's resource id (or --bucket when given explicitly).",
    )
    parser.add_argument(
        "--version",
        type=str,
        default=DEFAULT_DEPENDENCIES_VERSION,
        dest="dep_version",
        help=f"Dependency-bundle version for the bucket layout (default: {DEFAULT_DEPENDENCIES_VERSION}). "
        "Workflows are registered as <wdl>-v<version> (e.g. -v1-7-0). Must match the cohort's dependencies_version "
        "in the run config.",
    )
    parser.add_argument(
        "--repo-root",
        type=str,
        default=None,
        dest="repo_root",
        help="Path to the biobank_workflows checkout holding the WDLs / dependency scripts. "
        "Auto-detected from the fedcat package and cwd when omitted.",
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--force", action="store_true", help="Replace existing workflows without prompting.")
    group.add_argument(
        "--skip-existing",
        action="store_true",
        dest="skip_existing",
        help="Leave existing workflows untouched.",
    )
    parser.add_argument("--dry-run", action="store_true", dest="dry_run", help="Log actions without executing.")
    return parser


def _setup_workflows_main(args: argparse.Namespace) -> None:
    """Run the `setup-workflows` subcommand."""
    from fedcat.platforms.gcs_cromwell.buckets import resolve_workspace_bucket
    from fedcat.platforms.gcs_cromwell.workflow_setup import WorkflowSetupError, setup_workflows

    biobank = BioBank(args.biobank)
    if args.bucket:
        bucket = args.bucket
        bucket_id = args.bucket_id or args.bucket
    else:
        try:
            chosen = resolve_workspace_bucket()
        except (RuntimeError, ValueError) as exc:
            logger.error(str(exc))
            sys.exit(1)
        bucket = chosen.gs_url
        bucket_id = args.bucket_id or chosen.resource_id

    try:
        setup_workflows(
            biobank=biobank,
            bucket=bucket,
            bucket_id=bucket_id,
            version=args.dep_version,
            repo_root=args.repo_root,
            force=args.force,
            skip_existing=args.skip_existing,
            dry_run=args.dry_run,
        )
    except NotImplementedError as exc:
        logger.error(str(exc))
        sys.exit(1)
    except WorkflowSetupError as exc:
        logger.error(str(exc))
        sys.exit(1)


def _add_cromwell_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Add the `cromwell` subcommand's arguments."""
    from fedcat.platforms.gcs_cromwell.cromwell_server import DEFAULT_PORT

    parser.add_argument(
        "action",
        choices=["setup", "start", "status", "stop"],
        help="setup: write cromwell.conf via 'wb cromwell generate-config' plus fedcat's settings; "
        "start: launch the server in the background and wait until healthy; "
        "status: check the server is healthy and on a cloud backend; stop: shut it down.",
    )
    parser.add_argument(
        "--dir",
        type=Path,
        default=Path("~/cromwell"),
        dest="cromwell_dir",
        help="Directory for cromwell.conf, its log, pid file and database (default: ~/cromwell).",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=DEFAULT_PORT,
        help=f"Server port (default: {DEFAULT_PORT}). fedcat discover reads CROMWELL_URL when it is not the default.",
    )
    parser.add_argument(
        "--bucket-id",
        type=str,
        default=None,
        dest="bucket_id",
        help="setup: Workbench bucket resource id to root workflow execution in. "
        "If omitted, you'll choose from the workspace's CONTROLLED buckets (via wb).",
    )
    parser.add_argument(
        "--concurrent-job-limit",
        type=int,
        default=2000,
        dest="concurrent_job_limit",
        help="setup: most tasks Cromwell runs at once (default: 2000). Keep within the project's Batch quota.",
    )
    parser.add_argument("--force", action="store_true", help="setup: overwrite an existing cromwell.conf.")
    parser.add_argument(
        "--memory",
        type=str,
        default="16g",
        help="start: JVM heap for the server (default: 16g); raise for phases with many thousands of tasks.",
    )
    parser.add_argument(
        "--cromwell",
        type=str,
        default="cromwell",
        help="start: the cromwell command, or a path to cromwell.jar (default: cromwell).",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=300,
        help="start/stop: seconds to wait for the server to come up or shut down (default: 300).",
    )
    return parser


def _cromwell_main(args: argparse.Namespace) -> None:
    """Run the `cromwell` subcommand."""
    from fedcat.platforms.gcs_cromwell import cromwell_server as cs

    directory = args.cromwell_dir.expanduser()
    try:
        match args.action:
            case "setup":
                bucket_id = args.bucket_id
                if not bucket_id:
                    from fedcat.platforms.gcs_cromwell.buckets import resolve_workspace_bucket

                    bucket_id = resolve_workspace_bucket().resource_id
                conf = cs.write_config(directory, bucket_id, args.concurrent_job_limit, force=args.force)
                logger.success(f"Wrote {conf}; start the server with 'fedcat cromwell start'")
            case "start":
                status = cs.start_server(
                    directory, port=args.port, memory=args.memory, cromwell=args.cromwell, timeout=args.timeout
                )
                healthy, lines = cs.describe(status, args.port)
                for line in lines:
                    (logger.success if healthy else logger.error)(line)
                if args.port != cs.DEFAULT_PORT:
                    logger.info(f"Before running fedcat discover: export CROMWELL_URL={cs.server_url(args.port)}")
                if not healthy:
                    sys.exit(1)
            case "status":
                healthy, lines = cs.describe(cs.server_status(args.port), args.port)
                for line in lines:
                    (logger.success if healthy else logger.error)(line)
                if not healthy:
                    sys.exit(1)
            case "stop":
                if cs.stop_server(directory, timeout=args.timeout):
                    logger.success("Cromwell stopped; running workflows resume when it restarts")
                else:
                    logger.info(f"No Cromwell server started from {directory} is running")
    except (cs.CromwellServerError, RuntimeError, ValueError) as exc:
        logger.error(str(exc))
        sys.exit(1)


def _build_parser() -> argparse.ArgumentParser:
    """Build the `fedcat` argument parser, with every subcommand.

    Also used by `tools/docs/gen_docs.py` and the tests.
    """
    top_parser = argparse.ArgumentParser(
        prog="fedcat",
        description=(
            f"fedcat {__version__} — unified AoU + UKB GWAS discovery and "
            "forward-selection federated conditioning runner."
        ),
    )
    top_parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="Show program version and exit.",
    )
    subparsers = top_parser.add_subparsers(dest="subcommand")
    subparsers.required = True

    # ── discover subcommand ────────────────────────────────────────────────────
    disc_parser = subparsers.add_parser(
        "discover",
        description=(
            "Run GWAS discovery and federated conditioning. "
            "Conditioning parameters (p-value threshold, r² threshold etc.) "
            "live in the config JSON."
        ),
        help="Run GWAS discovery and federated conditioning.",
    )
    _add_discover_args(disc_parser)

    # ── convert-annots subcommand ─────────────────────────────────────────────
    conv_parser = subparsers.add_parser(
        "convert-annots",
        description=(
            "Convert annotation TSV files to the optimised parquet schema used "
            "by the ACAT engine.  Run this once per annotation set and pass "
            "--schema-dir to downstream ACAT/conditioning steps."
        ),
        help="Convert annotation TSV files to optimised parquet schema.",
    )
    _add_convert_annots_args(conv_parser)

    # ── meta subcommand ───────────────────────────────────────────────────────
    meta_parser = subparsers.add_parser(
        "meta",
        description=(
            "Run IVW fixed-effects meta-analysis on per-cohort REGENIE output "
            "files.  Input files are defined via a cohort file (see --cohort-file)."
        ),
        help="Run IVW fixed-effects meta-analysis.",
    )
    _add_meta_args(meta_parser)

    # ── acat subcommand ───────────────────────────────────────────────────────
    acat_parser = subparsers.add_parser(
        "acat",
        description=(
            "Compute ACAT p-values from single-variant and optional burden "
            "meta-analysis outputs.  Requires either --schema-dir (parquet) "
            "or --masks-tsvs + --annotations-tsvs."
        ),
        help="Compute ACAT p-values.",
    )
    _add_acat_args(acat_parser)

    # ── subset-results subcommand ────────────────────────────────────────────
    subset_parser = subparsers.add_parser(
        "subset-results",
        description=(
            "Subset meta-analysis results for conditioning (step 0). "
            "Scans burden, ACAT, and singles outputs for significant hits "
            "and packages tarballs for downstream conditioning."
        ),
        help="Subset meta-analysis results for conditioning (step 0).",
    )
    _add_subset_results_args(subset_parser)

    # ── summarise-indep subcommand ────────────────────────────────────────────
    summarise_parser = subparsers.add_parser(
        "summarise-indep",
        description=(
            "Regenerate conditioning summaries (all_sig / all_top_sig / trace files) "
            "for a finished or old run, from its per-step summary_output tables and a "
            "JSON containing top-level 'passing_ids'.  Back-compatible with old fedcat/"
            "fedcond/looper versions; no cloud upload or state persistence is performed."
        ),
        help="Regenerate conditioning summaries for a finished run.",
    )
    _add_summarise_indep_args(summarise_parser)

    # ── collinearity subcommand ──────────────────────────────────────────────
    collinearity_parser = subparsers.add_parser(
        "collinearity",
        description=(
            "Check collinearity between selected (conditioned-on) and candidate "
            "variants from a plink2 pgen/pvar. Reports per-candidate VIF R² against "
            "the selected set and the max VIF R² across the augmented set."
        ),
        help="Check collinearity of candidate variants against a selected set.",
    )
    _add_collinearity_args(collinearity_parser)

    # ── setup-workflows subcommand ────────────────────────────────────────────
    setup_wf_parser = subparsers.add_parser(
        "setup-workflows",
        description=(
            "Stage a biobank's WDL workflows and their runtime dependency scripts "
            "into the current Workbench workspace. Uploads dependencies to the "
            "workspace bucket, bakes the per-workspace dependencies_path into each "
            "WDL, and registers them via 'wb workflow create'. Only 'aou' is "
            "currently implemented."
        ),
        help="Stage repo WDL workflows + dependencies into the workspace.",
    )
    _add_setup_workflows_args(setup_wf_parser)

    # ── cromwell subcommand ───────────────────────────────────────────────────
    cromwell_parser = subparsers.add_parser(
        "cromwell",
        description=(
            "Manage a Cromwell server for the local_cromwell backend: write its config "
            "(Google Batch, rooted in a workspace bucket), start it in the background, "
            "check its health, or stop it."
        ),
        help="Set up, start, check or stop a Cromwell server for local_cromwell.",
    )
    _add_cromwell_args(cromwell_parser)

    return top_parser


def cli() -> None:
    """Run the `fedcat` command: parse the arguments and run the subcommand."""
    args = _build_parser().parse_args()

    try:
        match args.subcommand:
            case "discover":
                # --force-job-retry is --reset-jobs 0 --force-reset
                if args.force_job_retry:
                    args.reset_jobs = 0
                    args.force_reset = True
                if sum((args.abort_all, args.report, args.list_jobs, args.dry_run)) > 1:
                    logger.error("--abort-all, --report, --list-jobs and --dry-run cannot be combined.")
                    sys.exit(1)
                asyncio.run(_discover_main(args))
            case "convert-annots":
                _convert_annots_main(args)
            case "meta":
                _meta_main(args)
            case "acat":
                _acat_main(args)
            case "subset-results":
                _subset_results_main(args)
            case "summarise-indep":
                _summarise_indep_main(args)
            case "collinearity":
                _collinearity_main(args)
            case "setup-workflows":
                _setup_workflows_main(args)
            case "cromwell":
                _cromwell_main(args)
    except KeyboardInterrupt:
        if args.subcommand == "discover":
            logger.warning("Paused (any launched jobs are still running). Run the same command again to resume.")
        else:
            logger.warning("Interrupted by user, exiting...")
        sys.exit(130)


if __name__ == "__main__":
    cli()
