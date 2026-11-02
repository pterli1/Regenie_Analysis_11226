"""Orchestration of a federated discovery run.

[`DiscoveryRun`][fedcat.core.engine.DiscoveryRun] runs each phase in turn, from
step 1 to conditioning, and saves its progress to resume from.
"""

from __future__ import annotations

import asyncio
import json
import math
import re
import shutil
import tarfile
import weakref
from collections.abc import Awaitable, Callable
from datetime import datetime
from enum import Enum
from itertools import accumulate
from pathlib import Path
from typing import Any

import polars as pl
from attrs import Factory, define, field
from loguru import logger

from fedcat.core.job import Job
from fedcat.core.logging import (
    COND_STEP_VAR,
    PHENO_VAR,
    STAGE_VAR,
    log_context,
    periodic_progress,
    write_text_atomic,
)
from fedcat.core.persistence import PersistentState
from fedcat.core.structures import (
    BIOBANK_PROFILES,
    PHASE_ORDERED,
    BioBank,
    CondType,
    DiscoveryConfig,
    JobState,
    MetaAnalysisConfig,
    MetaMode,
    Phase,
    RunEnv,
    build_exclude_covars,
    make_conv,
    pin_legacy_config,
    rebase_path,
)
from fedcat.core.tracking import progress_line
from fedcat.platforms.base import Platform
from fedcat.platforms.gcs_cromwell import gcs
from fedcat.stats.subset_results import _merge_masks_tsv

# Bump when the persisted shape of DiscoveryRun, Job or ConditioningStep changes
# incompatibly; state files with another version are refused on load.
STATE_VERSION: int = 5


# ---------------------------------------------------------------------------
# Config comparison
# ---------------------------------------------------------------------------


class CondOutcome(Enum):
    """How a phenotype's conditioning ended, or how a conditioning phase did overall."""

    COMPLETED = "completed"
    # Stopped early on purpose (--cond-one-step, cond_max_steps); resuming carries on
    INCOMPLETE = "incomplete"
    FAILED = "failed"


# Set at runtime rather than by the config file, so never a difference
_DERIVED_CONFIG_PATHS = re.compile(r"^cohorts\.[^.]+\.output_location$")
# Limit how far a run goes rather than what it computes, so a resume takes the current value
RUN_CONTROL_FIELDS: tuple[str, ...] = ("cond_max_steps",)
_MISSING = object()


def config_differences(current: DiscoveryConfig, persisted: DiscoveryConfig) -> list[tuple[str, Any, Any]]:
    """List the config values that differ between two configs.

    Nested objects and mappings are compared down to their leaves; lists are compared
    whole. Derived values and `RUN_CONTROL_FIELDS` are ignored.

    Args:
        current: Config given on the command line.
        persisted: Config saved in the run state.

    Returns:
        `(dotted path, persisted value, current value)` per difference.
    """
    conv = make_conv()
    diffs: list[tuple[str, Any, Any]] = []

    def _walk(old: Any, new: Any, path: str) -> None:
        """Record the differing leaves under `path`."""
        if isinstance(old, dict) and isinstance(new, dict):
            for key in sorted(set(old) | set(new), key=str):
                _walk(old.get(key, _MISSING), new.get(key, _MISSING), f"{path}.{key}" if path else str(key))
        elif old != new and not _DERIVED_CONFIG_PATHS.match(path) and path not in RUN_CONTROL_FIELDS:
            diffs.append((path, old, new))

    _walk(conv.unstructure(persisted), conv.unstructure(current), "")
    return diffs


def adopt_run_controls(config: DiscoveryConfig, current: DiscoveryConfig) -> None:
    """Copy the `RUN_CONTROL_FIELDS` values of the current config into a persisted one.

    Args:
        config: Persisted config; updated in place.
        current: Config given on the command line.
    """
    for name in RUN_CONTROL_FIELDS:
        setattr(config, name, getattr(current, name))


def compare_configs(current: DiscoveryConfig, persisted: DiscoveryConfig) -> list[str]:
    """Return the dotted paths of the config values that differ.

    See [`config_differences`][fedcat.core.engine.config_differences].

    Args:
        current: Config given on the command line.
        persisted: Config saved in the run state.

    Returns:
        Dotted paths of the differing values.
    """
    return [path for path, _, _ in config_differences(current, persisted)]


def format_config_differences(diffs: list[tuple[str, Any, Any]]) -> str:
    """Format config differences as one `path: persisted → current` line each.

    Args:
        diffs: Differences from [`config_differences`][fedcat.core.engine.config_differences].

    Returns:
        The lines, with values over 80 characters shortened.
    """

    def _show(value: Any) -> str:
        """Return a value's repr, shortened to 80 characters."""
        text = "(unset)" if value is _MISSING else repr(value)
        return text if len(text) <= 80 else f"{text[:77]}..."

    return "\n".join(f"  {path}: {_show(old)} → {_show(new)}" for path, old, new in diffs)


# ---------------------------------------------------------------------------
# Block-file helpers
# ---------------------------------------------------------------------------


def _normalise_chrom(chrom: str) -> str:
    """Strip whitespace and any `chr` prefix from a chromosome name."""
    c = str(chrom).strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    return c


def _read_blockfile(path: str | Path) -> list[tuple[str, int, int]]:
    """Parse a whitespace-separated `CHROM START END` block file, with an optional header.

    Lines starting with `#` and lines with fewer than three columns are skipped, and
    `chr` prefixes are removed.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)

    out: list[tuple[str, int, int]] = []
    with p.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            # Skip a header row
            if not out:
                try:
                    int(parts[1])
                    int(parts[2])
                except ValueError:
                    p0, p1, p2 = parts[0].lower(), parts[1].lower(), parts[2].lower()
                    if p0 in ("chrom", "chr", "chromosome") and "start" in p1 and ("end" in p2 or "stop" in p2):
                        continue
                    raise
            chrom = _normalise_chrom(parts[0])
            start = int(parts[1])
            end = int(parts[2])
            if end < start:
                raise ValueError(f"END < START in {p}: {start}, {end}")
            out.append((chrom, start, end))

    if not out:
        raise ValueError(f"No blocks found in {p}")
    return out


def _unsupported_interface(interface: str, bb: BioBank | None = None) -> str:
    """Return the error message for a `workflow_interface` a `match` does not handle."""
    where = f" (biobank {bb.value!r})" if bb is not None else ""
    return f"unsupported workflow_interface {interface!r}{where}"


def _step2_block_length(
    interface: str, is_singles: bool, coding_only: bool, size: int, max_length: int | None = None
) -> int:
    """Return the target step 2 block length in bp.

    Longer for smaller cohorts, scaled per workflow interface, with a floor and, for
    aggregates, a ceiling. Coding-only aggregates use the ceiling.
    """
    match interface:
        case "wdl_scatter_v1":
            singles_scale, agg_scale, agg_ceiling = 1.5e12, 1.5e11, 20_000_000
        case "dx_applet_v1":
            singles_scale, agg_scale, agg_ceiling = 3e12, 3e11, 30_000_000
        case _:
            raise ValueError(_unsupported_interface(interface))
    if is_singles:
        length = max(singles_scale / size, 6_000_000)
    else:
        ceiling = max_length or agg_ceiling
        length = ceiling if coding_only else max(min(agg_scale / size, ceiling), 600_000)
    return int(min(length, max_length) if max_length else length)


def custom_blocks_array(
    var_blocks: list[tuple[str, int, int]],
    target_block_length: int,
    prefix: str | None,
    no_pgs_chroms: list[str],
    chroms_to_run: list[str],
    fill_gaps: bool = True,
    regions: list[tuple[str, int, int]] | None = None,
    exclude_covars: list[str] | None = None,
    exclude_covars_dict: dict[str, list[str]] | None = None,
) -> list[list[str]]:
    """Group variant-density blocks into REGENIE step 2 blocks of about a target length.

    Each chromosome's blocks are split into near-equal groups. With `regions`, only
    blocks overlapping each region are grouped, and the groups are clipped to it.

    Args:
        var_blocks: `(chrom, start, end)` blocks of about 100k variants each.
        target_block_length: Target step 2 block length in bp.
        prefix: PGS covariate prefix (see
            [`build_exclude_covars`][fedcat.core.structures.build_exclude_covars]).
        no_pgs_chroms: Chromosomes without a PGS covariate.
        chroms_to_run: Chromosomes to keep.
        fill_gaps: Extend blocks to cover the gaps between them, and each chromosome's
            first and last block to its ends (whole-genome mode only).
        regions: `(chrom, start, end)` regions to restrict testing to, or None.
        exclude_covars: Covariates excluded on every chromosome.
        exclude_covars_dict: Further covariates excluded per chromosome.

    Returns:
        `[chrom, start, stop, exclude_covars]` per block, all strings, as the workflow
        and applet inputs take them.
    """
    blocks_per_group = max(target_block_length // 100_000, 1)
    chroms_set = set(chroms_to_run)

    def _excl(chrom: str) -> str:
        """Return the covariates excluded on a chromosome."""
        return build_exclude_covars(
            chrom=chrom,
            pgs_prefix=prefix,
            chromosomes_no_pgs=no_pgs_chroms,
            exclude_covars=exclude_covars or (),
            exclude_covars_dict=exclude_covars_dict,
        )

    if not regions:
        chr_list = [r[0] for r in var_blocks]
        start_list = [r[1] for r in var_blocks]
        stop_list = [r[2] for r in var_blocks]

        group_sizes = [0]
        for chrom in dict.fromkeys(chr_list):
            num_blocks = chr_list.count(chrom)
            num_groups = math.ceil(num_blocks / blocks_per_group)
            base_size = num_blocks // num_groups
            extras = num_blocks % num_groups
            sizes = [base_size + 1 if i < extras else base_size for i in range(num_groups)]
            group_sizes += sizes
            assert sum(sizes) == num_blocks

        block_indexes = list(accumulate(group_sizes))
        idx_start = block_indexes[:-1]
        idx_stop = [x - 1 for x in block_indexes[1:]]

        chroms = [chr_list[i] for i in idx_start]
        starts = [start_list[i] + 1 for i in idx_start]
        starts[0] -= 1
        stops = [stop_list[i] for i in idx_stop]
        exclcovars = [_excl(ch) for ch in chroms]

        if fill_gaps:
            i = 0
            while i < len(chroms):
                chrom = chroms[i]
                j = i
                while j < len(chroms) and chroms[j] == chrom:
                    j += 1
                starts[i] = 1
                for k in range(i + 1, j):
                    starts[k] = min(starts[k], stops[k - 1] + 1)
                stops[j - 1] = 300_000_000
                i = j

        return [
            [str(ch), str(st), str(sp), xc]
            for ch, st, sp, xc in zip(chroms, starts, stops, exclcovars)
            if ch in chroms_set
        ]

    # Regions mode: blocks grouped within each region
    out_rows: list[list[str]] = []
    for chrom, region_start, region_end in regions:
        chrom = str(chrom)
        if chrom not in chroms_set:
            continue
        vb_chr = [
            (ch, st, sp)
            for ch, st, sp in var_blocks
            if ch == chrom and sp >= region_start and st + 1 <= region_end
        ]
        num_blocks = len(vb_chr)
        if num_blocks == 0:
            continue
        num_groups = math.ceil(num_blocks / blocks_per_group)
        base_size = num_blocks // num_groups
        extras = num_blocks % num_groups
        sizes = [base_size + 1 if i < extras else base_size for i in range(num_groups)]
        assert sum(sizes) == num_blocks
        block_indexes = list(accumulate([0] + sizes))
        idx_start = block_indexes[:-1]
        idx_stop = [x - 1 for x in block_indexes[1:]]
        starts_r = [vb_chr[i][1] + 1 for i in idx_start]
        stops_r = [vb_chr[i][2] for i in idx_stop]
        starts_r[0] = region_start
        stops_r[-1] = region_end
        excls = [_excl(chrom)] * num_groups
        if fill_gaps:
            for k in range(1, len(starts_r)):
                starts_r[k] = min(starts_r[k], stops_r[k - 1] + 1)
        out_rows.extend([chrom, str(s), str(e), x] for s, e, x in zip(starts_r, stops_r, excls))
    return out_rows


def _merge_step0_subsets(pheno: str, pheno_dl: Path, base_path: Path) -> None:
    """Merge one phenotype's per-chromosome step-0 subset tarballs.

    Writes `<pheno>_agg_cond_step0_subset.tar.gz` and
    `<pheno>_singles_cond_step0_subset.tar.gz` into `base_path`, for the types with
    tarballs in `pheno_dl`. Aggregate annotations may be flat TSVs or an
    `annotations_parquet/` schema directory.

    Raises:
        ValueError: If aggregate tophits exist without their annotation files.
    """
    from fedcat.stats.subset_results import _PARQUET_SCHEMA_FILES

    tarballs = sorted(pheno_dl.glob("*.tar.gz"))
    if not tarballs:
        logger.info("No step0 subset tarballs")
        return

    agg_dir = pheno_dl / f"{pheno}_agg_cond_step0_subset"
    sng_dir = pheno_dl / f"{pheno}_singles_cond_step0_subset"
    for d in (agg_dir, sng_dir):
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True, exist_ok=True)

    extracted = pheno_dl / "temp_extracted"
    if extracted.exists():
        shutil.rmtree(extracted)
    extracted.mkdir(parents=True, exist_ok=True)

    agg_tophits_paths: list[str] = []
    agg_annotations_paths: list[str] = []
    agg_masks_paths: list[Path] = []
    # Extracted paths per parquet schema file
    agg_parquet_paths: dict[str, list[str]] = {f: [] for f in _PARQUET_SCHEMA_FILES}
    sng_tophits_paths: list[str] = []
    has_parquet = False

    logger.debug(f"Extracting members from {len(tarballs)} tarballs")
    for tb in tarballs:
        with tarfile.open(tb, "r:gz") as tar_in:
            for member in tar_in.getmembers():
                name = Path(member.name).name
                parent = Path(member.name).parent.name if "/" in member.name else ""

                if parent == "annotations_parquet" and name in _PARQUET_SCHEMA_FILES:
                    out_path = extracted / f"{tb.stem}__{name}"
                    f = tar_in.extractfile(member)
                    if f is None:
                        logger.warning(f"Failed to extract {member.name} from {tb}")
                        continue
                    with open(out_path, "wb") as fh:
                        shutil.copyfileobj(f, fh)
                    agg_parquet_paths[name].append(str(out_path))
                    has_parquet = True
                    continue

                if member.isdir():
                    continue

                if name not in ("tophits_input", "annotations", "masks"):
                    logger.warning(f"Unexpected member '{name}' in {tb}, skipping")
                    continue
                out_path = extracted / f"{tb.stem}__{name}"
                f = tar_in.extractfile(member)
                if f is None:
                    logger.warning(f"Failed to extract {name} from {tb}, skipping")
                    continue
                with open(out_path, "wb") as fh:
                    shutil.copyfileobj(f, fh)
                if name == "tophits_input":
                    if "_agg_cond_step0_subset" in tb.name:
                        agg_tophits_paths.append(str(out_path))
                    elif "_singles_cond_step0_subset" in tb.name:
                        sng_tophits_paths.append(str(out_path))
                elif name == "annotations":
                    agg_annotations_paths.append(str(out_path))
                elif name == "masks":
                    agg_masks_paths.append(out_path)

    def _merge_paths_to_file(in_paths: list[str], out_path: Path) -> None:
        """Write the unique rows of headerless TSVs to `out_path`, or an empty file if none."""
        if not in_paths:
            logger.warning(f"No input files to merge for {out_path}, creating empty file")
            out_path.write_text("")
            return
        pl.scan_csv(
            in_paths,
            has_header=False,
            separator="\t",
            infer_schema=False,
            raise_if_empty=False,
        ).unique().sink_csv(out_path, include_header=False, separator="\t", engine="streaming")

    logger.debug("Merging extracted files")
    _merge_paths_to_file(agg_tophits_paths, agg_dir / "tophits_input")
    _merge_paths_to_file(sng_tophits_paths, sng_dir / "tophits_input")

    # Aggregate tophits need their annotation files
    if agg_tophits_paths:
        if has_parquet:
            missing_schema = [f for f, paths in agg_parquet_paths.items() if not paths]
            if missing_schema:
                raise ValueError(
                    "Missing one or more annotations_parquet schema files while aggregate tophits exist "
                    f"for {pheno}: {', '.join(sorted(missing_schema))}"
                )
        else:
            missing_flat = []
            any_tophits = len((agg_dir / "tophits_input").read_text().strip()) > 0
            if not agg_annotations_paths and any_tophits:
                missing_flat.append("annotations")
            if not agg_masks_paths and any_tophits:
                missing_flat.append("masks")
            if missing_flat:
                raise ValueError(
                    "Missing flat annotation files while aggregate tophits exist "
                    f"for {pheno}: {', '.join(missing_flat)}"
                )

    if has_parquet:
        merged_pq_dir = agg_dir / "annotations_parquet"
        merged_pq_dir.mkdir(parents=True, exist_ok=True)
        for schema_file, paths in agg_parquet_paths.items():
            merged_path = merged_pq_dir / schema_file
            pl.concat([pl.read_parquet(p) for p in paths]).unique().write_parquet(merged_path)
        logger.debug("Merged parquet schema")
    else:
        _merge_paths_to_file(agg_annotations_paths, agg_dir / "annotations")
        _merge_masks_tsv(agg_masks_paths, agg_dir / "masks")

    logger.debug("Creating final subset tarballs")
    if agg_tophits_paths or agg_annotations_paths or has_parquet:
        with tarfile.open(base_path / f"{pheno}_agg_cond_step0_subset.tar.gz", "w:gz") as tar_out:
            tar_out.add(agg_dir / "tophits_input", arcname="tophits_input")
            if has_parquet:
                tar_out.add(agg_dir / "annotations_parquet", arcname="annotations_parquet")
            else:
                tar_out.add(agg_dir / "annotations", arcname="annotations")
                tar_out.add(agg_dir / "masks", arcname="masks")
    else:
        logger.info("No aggregate subset parts; skipping agg subset tarball")

    if sng_tophits_paths:
        with tarfile.open(base_path / f"{pheno}_singles_cond_step0_subset.tar.gz", "w:gz") as tar_out:
            tar_out.add(sng_dir / "tophits_input", arcname="tophits_input")
    else:
        logger.info("No unconditional singles subset parts; skipping singles subset tarball")

    logger.info("Merged conditioning subset tarballs")


# ---------------------------------------------------------------------------
# DiscoveryRun
# ---------------------------------------------------------------------------


@define
class DiscoveryRun(PersistentState):
    """A discovery run: runs the requested phases in order, saving state to resume from.

    Each phase builds its jobs and runs them; the conditioning phases run a
    [`ConditioningStep`][fedcat.core.conditioning.ConditioningStep] loop per phenotype.
    State is saved per phase to `discovery_state_<phase>.json` in `run_path`.
    """

    config: DiscoveryConfig
    """Run config."""
    run_path: Path
    """Local run directory, holding state files, logs and downloads."""
    variant_density_blockfile: Path
    """Block file of about 100k variants per block, used to size step 2 blocks."""
    state_version: int = STATE_VERSION
    """Persisted state version, checked on load."""
    meta_mode: MetaMode = MetaMode.ALL
    """Result types to meta-analyse."""
    run_blockfile: Path | None = None
    """Regions to restrict step 2 to, or None for whole chromosomes."""
    phases_to_run: list[Phase] = field(factory=lambda: PHASE_ORDERED.copy())
    """Phases to run, in order."""
    current_phase: Phase = Phase.STEP1
    """Phase being run."""
    required: bool = True
    """Whether phases remain to run."""
    had_error: bool = False
    """Whether the current phase stopped with an error."""
    jobs: list[Job] = Factory(list)
    """Jobs of the current phase."""
    cond_status: dict[CondType, dict[str, CondOutcome]] = field(factory=dict)
    """Outcome of each phenotype's conditioning, per conditioning type."""

    platforms: dict[BioBank, Platform] = field(init=False, factory=dict, repr=False)
    """Platform per biobank, rebuilt from the config on load; not persisted."""

    # Runtime controls set from the command line; not persisted
    cond_step: int | None = field(default=None, metadata={"cattrs_exclude": True})
    """Conditioning step to start from; None to detect it from the step directories."""
    cond_one_step: bool = field(default=False, metadata={"cattrs_exclude": True})
    """Run one conditioning step, then stop."""
    cond_reset_step: int | None = field(default=None, metadata={"cattrs_exclude": True})
    """`StepState` value to reset a failed conditioning step to."""
    cond_reset_jobs: int | None = field(default=None, metadata={"cattrs_exclude": True})
    """`JobState` value to reset failed conditioning jobs to."""
    cond_force_reset: bool = field(default=False, metadata={"cattrs_exclude": True})
    """Force a conditioning step into ERROR before applying `cond_reset_step`."""
    cond_retry_failed: bool = field(default=False, metadata={"cattrs_exclude": True})
    """Resubmit failed conditioning jobs, or rerun the stage a step failed in (`--force-job-retry`)."""
    cond_singles_tarball: Path | None = field(default=None, metadata={"cattrs_exclude": True})
    """Step-0 singles subset tarball to use instead of the one in the run directory."""
    cond_agg_tarball: Path | None = field(default=None, metadata={"cattrs_exclude": True})
    """Step-0 aggregate subset tarball to use instead of the one in the run directory."""
    _overwrite_config: bool = field(default=False, metadata={"cattrs_exclude": True})
    """Replace the config saved in conditioning state with the current one (`--overwrite-config`)."""
    _stop_requested: bool = field(default=False, metadata={"cattrs_exclude": True})
    """Stop after the current phase without advancing (set by `--cond-one-step`)."""

    def __attrs_post_init__(self) -> None:
        """Create the run directory, link jobs to this run, and build the platforms."""
        from fedcat.platforms import build_platforms

        self.run_path.mkdir(parents=True, exist_ok=True)
        for job in self.jobs:
            job._parent = weakref.ref(self)

        for biobank in self.config.cohorts:
            self.config.cohorts[biobank].output_location = self.config.output_location(biobank)

        self.platforms = build_platforms(self.config)

    def mark_jobs_stale(self) -> None:
        """Mark every job's inputs to be rebuilt from the config before it is submitted."""
        for job in self.jobs:
            job._input_stale = True

    def adopt_config(self, config: DiscoveryConfig) -> None:
        """Replace the persisted config, and everything derived from it (`--overwrite-config`).

        `run_path` is re-derived and files under the old `local_run_path` are moved to
        the new one, so a moved run directory resumes in place. Job inputs are rebuilt
        from the new config before submission.

        Args:
            config: New config.
        """
        old_root, old_run_path = Path(self.config.local_run_path), self.run_path
        self.config = config
        self.run_path = config.discovery_run_path
        self.run_blockfile = rebase_path(self.run_blockfile, old_root, Path(config.local_run_path))
        self.variant_density_blockfile = rebase_path(  # type: ignore[assignment]
            self.variant_density_blockfile, old_root, Path(config.local_run_path)
        )
        if self.run_path != old_run_path:
            logger.info(f"Run directory moved with the config: {old_run_path} → {self.run_path}")
        self.__attrs_post_init__()
        self.mark_jobs_stale()

    # ------------------------------------------------------------------
    # Phase navigation
    # ------------------------------------------------------------------

    def advance_phase(self) -> None:
        """Move to the next phase in `phases_to_run`, clearing the jobs, or finish after the last."""
        try:
            idx = self.phases_to_run.index(self.current_phase)
        except ValueError:
            logger.error(f"Current phase {self.current_phase!r} not in phases_to_run {self.phases_to_run}")
            self.required = False
            self.had_error = True
            return
        if idx < len(self.phases_to_run) - 1:
            self.current_phase = self.phases_to_run[idx + 1]
            self.jobs.clear()
        else:
            self.required = False

    async def _job_progress_msg(self) -> str:
        """Return the progress line of the current jobs, for the periodic progress log."""
        return await progress_line(self.jobs, self.platforms)

    # ------------------------------------------------------------------
    # GCS / cloud helpers
    # ------------------------------------------------------------------

    def _meta_base_folder_gcs(self) -> str:
        """Return the GCS meta-analysis folder (without `gs://`), or `""` if there is none.

        Also locates the uploaded state files, so without a `meta_analysis_config` it
        follows `local_env`. A UKB-only run driven from AoU (with `aou_bucket` set) keeps
        them in the bucket, under the UKB cohort's output path.
        """
        meta_cfg = self.config.meta_analysis_config
        location = meta_cfg.location if meta_cfg is not None else self.config.local_env
        if location == RunEnv.AOU and BioBank.AOU in self.config.cohorts:
            return f"{self.config.output_location(BioBank.AOU)}/meta_analysis"
        if self.config.aou_bucket and BioBank.UKB in self.config.cohorts:
            bucket = self.config.aou_bucket.removeprefix("gs://")
            return f"{bucket}{self.config.output_location(BioBank.UKB)}/meta_analysis"
        return ""

    def _state_cloud_folder(self) -> str:
        """Return the GCS folder state files are uploaded to, or `""` for none."""
        base = self._meta_base_folder_gcs()
        if not base:
            return ""
        return f"{base}/orchestrator_state/discovery"

    # ------------------------------------------------------------------
    # Serialisation
    # ------------------------------------------------------------------

    def _state_file(self) -> Path:
        """Return the state file of the current phase."""
        return self.run_path / f"discovery_state_{self.current_phase.value}.json"

    def _state_queue_key(self) -> str:
        """Return the key this run's state writes are queued under."""
        return f"discovery_{self.config.run_prefix}"

    def _state_upload_destination(self, path: Path) -> str | None:
        """Return the GCS path to upload a state file to, or None to keep it local."""
        cloud_folder = self._state_cloud_folder()
        return f"{cloud_folder}/{path.name}" if cloud_folder else None

    async def _snapshot_state(self, reason: str) -> None:
        """Write a timestamped copy of the state to `snapshots/`, and upload it if possible."""
        conv = make_conv()
        ts = datetime.now().strftime("%Y%m%dT%H%M%SZ")
        path = self.run_path / "snapshots" / f"snapshot_{self.current_phase.value}_{reason}_{ts}.json"
        write_text_atomic(path, conv.dumps(self, indent=2))
        logger.info(f"Snapshot saved: {path}")
        cloud_folder = self._state_cloud_folder()
        if cloud_folder:
            destination = f"{cloud_folder}/{path.name}"
            try:
                await gcs.upload_file_async(path, destination=destination)
            except Exception as exc:  # noqa: BLE001 — a snapshot that fails to upload is only warned about
                logger.warning(f"Failed to upload snapshot {path.name} to gs://{destination}: {exc}")

    # ------------------------------------------------------------------
    # Main run loop
    # ------------------------------------------------------------------

    async def run(
        self,
        one_phase: bool = False,
        reset_jobs: JobState | None = None,
        reset_jobs_list: list[str] | None = None,
        force_reset: bool = False,
    ) -> None:
        """Run phases from `current_phase` until all are done or one fails.

        Jobs are only reset when the previous run stopped with an error, unless
        `force_reset` is given.

        Args:
            one_phase: Stop after the current phase.
            reset_jobs: State to reset failed jobs to; None returns them to their state
                before the error.
            reset_jobs_list: Descriptions of the jobs to reset; None resets every failed job.
            force_reset: Run, and reset jobs, even if the previous run did not stop with
                an error.
        """
        if force_reset:
            self.required = True
        attempt_job_reset = False
        if self.had_error:
            logger.warning(
                f"The previous run stopped with an error in phase {self.current_phase.value}; retrying"
            )
            self.had_error = False
            attempt_job_reset = True
        elif reset_jobs is not None or reset_jobs_list is not None:
            if not force_reset:
                logger.error(
                    "--reset-jobs(-list) given but the previous run did not stop with an error; "
                    "use --force-reset to override"
                )
                return
            attempt_job_reset = True

        while self.required and not self.had_error:
            STAGE_VAR.set(self.current_phase.label)
            try:
                await self._run_phase(
                    reset_jobs=reset_jobs,
                    reset_jobs_list=reset_jobs_list or [],
                    attempt_job_reset=attempt_job_reset,
                )
                if self._stop_requested:
                    logger.info("Stopping before the next phase as requested; resume to carry on")
                    break
                if not self.had_error:
                    logger.success("Phase completed successfully")
                    if one_phase:
                        break
                    self.persist_state(reason=f"phase_{self.current_phase.value}_complete", max_delay=0)
                    await self.flush_state()
                    self.advance_phase()
            except asyncio.CancelledError:
                logger.info("Saving state before pausing…")
                raise
            except Exception:
                logger.exception("Unexpected error while running phase")
                self.had_error = True
            finally:
                self.persist_state()
                await self.flush_state()

        if not self.required and not self.had_error:
            logger.success("All requested phases completed successfully")

    # ------------------------------------------------------------------
    # Phase dispatch
    # ------------------------------------------------------------------

    async def _run_phase(
        self,
        reset_jobs: JobState | None,
        reset_jobs_list: list[str],
        attempt_job_reset: bool = False,
    ) -> None:
        """Run the current phase, setting `had_error` if it fails."""
        logger.info("Starting phase")
        resets = (reset_jobs, reset_jobs_list, attempt_job_reset)

        # ── REGENIE PHASES ────────────────────────────────────────────
        if self.current_phase in (Phase.STEP1, Phase.STEP2_SINGLES, Phase.STEP2_AGG):
            await self._run_jobs(self._build_jobs_for_phase, *resets)
            if (
                not self.had_error
                and self.current_phase in (Phase.STEP2_SINGLES, Phase.STEP2_AGG)
                and not await self._verify_step2_output_lists(
                    is_singles=self.current_phase == Phase.STEP2_SINGLES
                )
            ):
                logger.error("STEP2 output list verification failed")
                self.had_error = True

        # ── META-ANALYSIS ─────────────────────────────────────────────
        elif self.current_phase == Phase.META:
            meta_cfg = self.config.meta_analysis_config
            if meta_cfg is None:
                logger.error("meta_analysis_config must be set in the config to run the meta-analysis phase")
                self.had_error = True
                return

            # Staging runs once, before building the jobs; a resume reuses the jobs as they were
            async def _stage_and_build() -> list[Job]:
                """Stage other cohorts' results for an AoU meta-analysis, then build its jobs."""
                if meta_cfg.location == RunEnv.AOU:
                    logger.info("Staging remote-cohort results into the meta bucket")
                    await self._stage_cohort_results_to_meta()
                return await self._build_jobs_for_phase()

            await self._run_jobs(_stage_and_build, *resets)

        # ── STEP-0 SUBSET ─────────────────────────────────────────────
        elif self.current_phase == Phase.STEP0_SUBSET:
            await self._run_step0_subset(*resets)

        # ── CONDITIONING PHASES ───────────────────────────────────────
        elif self.current_phase in (Phase.COND_SINGLES, Phase.COND_AGG):
            cond_type = CondType.SINGLEVAR if self.current_phase == Phase.COND_SINGLES else CondType.AGGREGATE
            outcome = await self._run_conditioning(cond_type)
            if outcome is CondOutcome.FAILED:
                self.had_error = True
            elif outcome is CondOutcome.INCOMPLETE and self.cond_one_step:
                self._stop_requested = True
            # The command-line conditioning controls apply to the first conditioning phase only
            self.cond_step = None
            self.cond_reset_step = None
            self.cond_reset_jobs = None
            self.cond_force_reset = False
            self.cond_retry_failed = False
        else:
            logger.error("Phase not yet implemented")
            self.had_error = True

    async def _run_jobs(
        self,
        build: Callable[[], Awaitable[list[Job]]],
        reset_jobs: JobState | None,
        reset_jobs_list: list[str],
        attempt_job_reset: bool,
    ) -> None:
        """Build the phase's jobs, or reset and resume those already built, then run them.

        Args:
            build: Builds the phase's jobs; called when there are none.
            reset_jobs: State to reset failed jobs to; None for their state before the error.
            reset_jobs_list: Descriptions of the jobs to reset; empty to reset every failed job.
            attempt_job_reset: Whether to reset jobs at all.
        """
        if not self.jobs:
            self.jobs = await build()
            self.persist_state(reason="jobs_built", max_delay=0)
        elif reset_jobs_list and attempt_job_reset:
            matched: list[str] = []
            for j in self.jobs:
                if j.description in reset_jobs_list:
                    if j.state in (JobState.ERROR, JobState.EMPTY, JobState.READY):
                        j.reset_error_state(reset_jobs)
                        j.id = ""
                        matched.append(j.description)
                    elif j.state == JobState.RUNNING:
                        logger.warning(f"Skipping reset of RUNNING job {j.description!r}")
                    elif j.state == JobState.COMPLETED:
                        logger.warning(f"Skipping reset of COMPLETED job {j.description!r}")
            if len(matched) != len(reset_jobs_list):
                unmatched = set(reset_jobs_list) - set(matched)
                logger.warning(f"--reset-jobs-list: no match for {unmatched}")
            logger.info(f"Reset {len(matched)} jobs in phase {self.current_phase.value!r}")
            self.persist_state(reason="jobs_reset_list", max_delay=0)
            if any(j.state == JobState.ERROR for j in self.jobs):
                errored = [j.description for j in self.jobs if j.state == JobState.ERROR]
                logger.warning(f"Jobs still in ERROR: {errored}")
        elif attempt_job_reset:
            failed = [j for j in self.jobs if j.state == JobState.ERROR]
            for j in failed:
                j.reset_error_state(reset_jobs)
            if failed:
                logger.info(f"Retrying {len(failed)} failed job(s)")
            else:
                logger.info("No failed jobs to retry; resuming")
            self.persist_state(reason="jobs_reset", max_delay=0)
        else:
            logger.info("Resuming from existing job states")

        await self._skip_completed_units()
        # Running jobs first, so they take the concurrency slots they already occupy.
        runnable_jobs = sorted(
            (j for j in self.jobs if j.state in (JobState.EMPTY, JobState.READY, JobState.RUNNING)),
            key=lambda j: j.state != JobState.RUNNING,
        )
        if runnable_jobs:
            logger.info(f"Running {len(runnable_jobs)}/{len(self.jobs)} job(s)")
        else:
            logger.info(f"No runnable jobs; all {len(self.jobs)} already completed or errored")
        slots = {
            bb: asyncio.Semaphore(cohort.max_concurrent_jobs)
            for bb, cohort in self.config.cohorts.items()
            if cohort.max_concurrent_jobs > 0
        }

        async def _run_job(job: Job, delay: float) -> None:
            """Run a job after `delay` seconds, within its biobank's concurrent job limit."""
            if (slot := slots.get(job.biobank)) is None:
                return await job.run(delay_seconds=delay)
            await asyncio.sleep(delay)
            async with slot:
                await job.run()

        tasks = [asyncio.create_task(_run_job(j, float(i))) for i, j in enumerate(runnable_jobs)]
        if tasks:
            async with periodic_progress(self._job_progress_msg):
                await asyncio.gather(*tasks)
        if any(j.state == JobState.ERROR for j in self.jobs):
            logger.error("Some jobs failed")
            self.had_error = True

    # ------------------------------------------------------------------
    # Job construction
    # ------------------------------------------------------------------

    async def _build_jobs_for_phase(self) -> list[Job]:
        """Build the jobs of the current phase, with their inputs."""
        cfg = self.config
        jobs: list[Job] = []

        def _make_job(job: Job) -> Job:
            """Link a job to this run."""
            job._parent = weakref.ref(self)
            return job

        # ── STEP 1 ────────────────────────────────────────────────────
        if self.current_phase == Phase.STEP1:
            for bb, cohort in cfg.cohorts.items():
                for job in self._build_step1_jobs(bb, cohort):
                    _make_job(job)
                    await job.generate_input()
                    jobs.append(job)

        # ── STEP 2 ────────────────────────────────────────────────────
        elif self.current_phase in (Phase.STEP2_SINGLES, Phase.STEP2_AGG):
            is_singles = self.current_phase == Phase.STEP2_SINGLES
            var_blocks = _read_blockfile(self.variant_density_blockfile)
            regions = _read_blockfile(self.run_blockfile) if self.run_blockfile else None
            for bb, cohort in cfg.cohorts.items():
                for job in self._build_step2_jobs(bb, cohort, is_singles, var_blocks, regions):
                    _make_job(job)
                    await job.generate_input()
                    jobs.append(job)

        elif self.current_phase == Phase.META and cfg.meta_analysis_config is not None:
            jobs = await self._build_meta_job(cfg.meta_analysis_config)

        elif self.current_phase == Phase.STEP0_SUBSET:
            jobs = await self._build_step0_jobs()

        return jobs

    @staticmethod
    def _wdl_jobs(bb: BioBank, cohort, op: str, units: list[tuple], description: str) -> list[Job]:
        """Return WDL jobs over `units`, split into workflows of at most `max_tasks_per_workflow` tasks.

        Units are ordered by phenotype set, so each workflow carries few sets' inputs.
        """
        size = cohort.max_tasks_per_workflow or len(units) or 1
        chunks = [units[i : i + size] for i in range(0, len(units), size)]
        return [
            Job(
                biobank=bb,
                description=description if len(chunks) == 1 else f"{description}_part{i + 1}of{len(chunks)}",
                recipe={"op": op, "units": [list(u) for u in chunk]},
            )
            for i, chunk in enumerate(chunks)
        ]

    def _build_step1_jobs(self, bb: BioBank, cohort) -> list[Job]:
        """Return a cohort's step 1 jobs: WDL workflows over all units, or one applet job per set and stratum."""
        cfg = self.config
        interface = BIOBANK_PROFILES[bb].workflow_interface
        n_sets = len(cfg.step1_pheno_sets)
        units = [(s, stratum) for s in range(n_sets) for stratum in cohort.strata]
        match interface:
            case "wdl_scatter_v1":
                return self._wdl_jobs(
                    bb, cohort, "step1", units, f"{cfg.run_prefix}_{bb.value}_regenie_step1_all_strata"
                )
            case "dx_applet_v1":
                return [
                    Job(
                        biobank=bb,
                        description=(
                            f"{cfg.run_prefix}_{bb.value}_{stratum}{f'_set{s + 1}' if n_sets > 1 else ''}_regenie_step1"
                        ),
                        recipe={"op": "step1", "stratum": stratum, "set": s},
                        priority="high",
                    )
                    for s, stratum in units
                ]
            case _:
                raise ValueError(_unsupported_interface(interface, bb))

    def _build_step2_jobs(self, bb: BioBank, cohort, is_singles: bool, var_blocks, regions) -> list[Job]:
        """Return a cohort's step 2 jobs: WDL workflows over all units, or one applet job per group, stratum and block.

        Block coordinates are stored in the recipes, so a resume keeps the blocks the
        jobs were built with.
        """
        cfg = self.config
        interface = BIOBANK_PROFILES[bb].workflow_interface
        op = "step2_singles" if is_singles else "step2_agg"
        max_length = cfg.max_step2_singles_block_length if is_singles else cfg.max_step2_aggregates_block_length

        def _blocks(stratum: str) -> list:
            """Step 2 blocks for a stratum, sized by its cohort size."""
            return custom_blocks_array(
                var_blocks=var_blocks,
                target_block_length=_step2_block_length(
                    interface, is_singles, cfg.coding_only, cohort.cohort_sizes[stratum], max_length
                ),
                prefix=cfg.pgs_prefix,
                no_pgs_chroms=cfg.chromosomes_no_pgs,
                chroms_to_run=cfg.chromosomes_to_run,
                fill_gaps=True,
                regions=regions,
                exclude_covars=cfg.exclude_covars,
                exclude_covars_dict=cfg.exclude_covars_dict,
            )

        kind = "singles" if is_singles else "aggregates"
        blocks = {stratum: _blocks(stratum) for stratum in cohort.strata}
        match interface:
            case "wdl_scatter_v1":
                units = [
                    (g, stratum, chrom, start, stop)
                    for g in range(len(cfg.step2_pheno_groups))
                    for stratum in cohort.strata
                    for chrom, start, stop, _ in blocks[stratum]
                ]
                return self._wdl_jobs(
                    bb, cohort, op, units, f"{cfg.run_prefix}_{bb.value}_{kind}_all_blocks_all_strata"
                )
            case "dx_applet_v1":
                jobs: list[Job] = []
                for g in range(len(cfg.step2_pheno_groups)):
                    label = cfg.step2_group_label(g)
                    for stratum in cohort.strata:
                        for chrom, start, stop, excl in blocks[stratum]:
                            jobs.append(
                                Job(
                                    biobank=bb,
                                    description=(
                                        f"{cfg.run_prefix}_{bb.value}_{stratum}{f'_{label}' if label else ''}"
                                        f"_{kind}_chr{chrom}_{start}_{stop}"
                                    ),
                                    recipe={
                                        "op": op,
                                        "stratum": stratum,
                                        "block": [chrom, start, stop, excl],
                                        "group": g,
                                    },
                                )
                            )
                return jobs
            case _:
                raise ValueError(_unsupported_interface(interface, bb))

    # ------------------------------------------------------------------
    # Staging results for the meta-analysis
    # ------------------------------------------------------------------

    def _meta_site_biobank(self) -> BioBank | None:
        """Return the first cohort on a `gcs_cromwell` platform, whose bucket a WDL meta-analysis reads."""
        for bb in self.config.cohorts:
            if BIOBANK_PROFILES[bb].platform_kind == "gcs_cromwell":
                return bb
        return None

    async def _stage_cohort_results_to_meta(self) -> bool:
        """Copy other cohorts' step 2 results into the meta-analysis bucket.

        Each cohort's results are downloaded from its platform and uploaded to
        `<meta_base>/<biobank>_results/<subdir>`, where the meta-analysis WDL reads them.
        """
        meta_site = self._meta_site_biobank()
        if meta_site is None:
            logger.info("No gcs_cromwell meta site; nothing to stage")
            return True
        meta_platform = self.platforms[meta_site]
        meta_base = self._meta_base_folder_gcs()
        subdir_names = ["single_variants"]
        if self.meta_mode.runs_aggregates:
            subdir_names.append("aggregates")
        pheno_pattern = "|".join(self.config.pheno_names)
        for bb, cohort in self.config.cohorts.items():
            if bb == meta_site:
                continue
            src_platform = self.platforms[bb]
            for subdir in subdir_names:
                stage_dir = self.run_path / ".meta_stage" / bb.value / subdir
                stage_dir.mkdir(parents=True, exist_ok=True)
                logger.info(f"Staging {bb.value} {subdir} results to {stage_dir}")
                await src_platform.download(
                    folder=f"{cohort.output_location}/meta_analysis/{subdir}",
                    pattern=rf".*_({pheno_pattern})\.regenie\.gz$",
                    dest=stage_dir,
                    max_workers=4,
                )
                dst = f"{meta_base}/{bb.value}_results/{subdir}"
                files = list(stage_dir.glob("*.regenie.gz"))
                logger.info(f"Uploading {len(files)} {bb.value} {subdir} files to gs://{dst}")
                for f in files:
                    await meta_platform.upload(f, destination=dst)
        return True

    def _step2_chunk_prefix(self, bb: BioBank, is_singles: bool, unit: tuple) -> str:
        """Return the name stem of one step 2 task's outputs, as the step 2 scripts build it."""
        group, anc, chrom, start, stop = unit
        label = self.config.step2_group_label(group)
        base = f"{self.config.run_prefix}_{anc}_chr{chrom}_{start}_{stop}{f'_{label}' if label else ''}"
        if is_singles:
            return f"{base}_step2_singlevar"
        annots = "_".join(e.name for e in self.config.cohorts[bb].regenie_annotations)
        return f"{base}_{annots}_step2_aggregates"

    def _expected_step2_output_lists(self, is_singles: bool) -> dict[BioBank, set[str]]:
        """Return the `.output_filelist` names each biobank's step 2 tasks should write."""
        op = "step2_singles" if is_singles else "step2_agg"
        expected: dict[BioBank, set[str]] = {bb: set() for bb in self.config.cohorts}
        for job in self.jobs:
            if job.recipe.get("op") == op:
                expected[job.biobank].update(
                    f"{self._step2_chunk_prefix(job.biobank, is_singles, u)}.output_filelist" for u in job.units
                )
        return expected

    def _completion_markers(self, bb: BioBank) -> tuple[str, str]:
        """Return the folder and filename pattern of the logs this phase's tasks upload last on success."""
        cohort = self.config.cohorts[bb]
        match self.current_phase:
            case Phase.STEP1:
                return self.platforms[bb].results_root(cohort.step1_directory), r".*_Step1_.*\.workflow\.log$"
            case Phase.STEP2_SINGLES:
                return f"{cohort.output_location}/meta_analysis/single_variants", r".*\.workflow\.log$"
            case _:
                return f"{cohort.output_location}/meta_analysis/aggregates", r".*\.workflow\.log$"

    def _unit_marker(self, job: Job, unit: tuple) -> str:
        """Return the completion log name of one step 1 or step 2 task."""
        cfg = self.config
        if job.recipe.get("op") == "step1":
            s, anc = unit
            prefix = cfg.step1_run_prefixes[s]
            return f"{prefix}_{cfg.trait_type.value}_{cfg.trait_transform.value}_Step1_{anc}.workflow.log"
        return f"{self._step2_chunk_prefix(job.biobank, job.recipe['op'] == 'step2_singles', unit)}.workflow.log"

    async def _skip_completed_units(self) -> None:
        """Remove already completed units from each unsubmitted job.

        Tasks also check their own completion log, but only once their machine has
        started; this avoids starting those machines. A job with nothing left is marked
        COMPLETED without being submitted.
        """
        pending = [j for j in self.jobs if j.state in (JobState.EMPTY, JobState.READY) and j.units]
        done: dict[BioBank, set[str]] = {}
        for bb in {j.biobank for j in pending}:
            folder, pattern = self._completion_markers(bb)
            try:
                done[bb] = set(await self.platforms[bb].list_files(folder=folder, pattern=pattern))
            except Exception as exc:  # noqa: BLE001 — platform errors are backend-specific; submitting everything is safe
                logger.warning(f"Could not list completed {bb.value} tasks in {folder}; submitting all: {exc}")
                done[bb] = set()
        for job in pending:
            remaining = [u for u in job.units if self._unit_marker(job, u) not in done[job.biobank]]
            if not remaining:
                logger.info(
                    f"All {len(job.units)} task(s) of {job.description!r} already complete; not submitting"
                )
                job.set_state(JobState.COMPLETED)
            elif remaining != job.pending_units:
                logger.info(
                    f"{job.description!r}: {len(job.units) - len(remaining)} of {len(job.units)} task(s) "
                    "already complete; submitting the rest"
                )
                job.recipe["pending"] = [list(u) for u in remaining]
                job._input_stale = True
        if pending:
            self.persist_state(reason="completed_units_skipped", max_delay=0)

    async def _verify_step2_output_lists(self, is_singles: bool) -> bool:
        """Check every step 2 task wrote its `.output_filelist`, and that the files it lists exist.

        Returns:
            True if all outputs are present; missing ones are logged.
        """
        subdir = "single_variants" if is_singles else "aggregates"
        expected_by_biobank = self._expected_step2_output_lists(is_singles=is_singles)
        check_root = self.run_path / ".step2_output_list_checks" / self.current_phase.value
        check_root.mkdir(parents=True, exist_ok=True)

        missing_lists: dict[BioBank, list[str]] = {}
        missing_outputs: dict[str, list[str]] = {}

        for biobank, expected_lists in expected_by_biobank.items():
            if not expected_lists:
                continue
            folder = f"{self.config.cohorts[biobank].output_location}/meta_analysis/{subdir}"
            local_dir = check_root / biobank.value
            local_dir.mkdir(parents=True, exist_ok=True)

            platform = self.platforms[biobank]
            remote_files = set(await platform.list_files(folder=folder))
            await platform.download(
                folder=folder,
                pattern=r".*\.output_filelist$",
                dest=local_dir,
                max_workers=4,
            )

            remote_lists = {name for name in remote_files if name.endswith(".output_filelist")}
            missing_expected_lists = sorted(expected_lists - remote_lists)
            if missing_expected_lists:
                missing_lists[biobank] = missing_expected_lists

            for manifest in sorted(expected_lists & remote_lists):
                manifest_path = local_dir / manifest
                if not manifest_path.exists():
                    missing_outputs[f"{biobank.value}:{manifest}"] = ["manifest_download_missing"]
                    continue
                listed_outputs = [line.strip() for line in manifest_path.read_text().splitlines() if line.strip()]
                missing_manifest_outputs = sorted({name for name in listed_outputs if name not in remote_files})
                if missing_manifest_outputs:
                    missing_outputs[f"{biobank.value}:{manifest}"] = missing_manifest_outputs

        if missing_lists or missing_outputs:
            for biobank, manifests in missing_lists.items():
                logger.error(
                    f"Missing {len(manifests)} expected .output_filelist files for {biobank.value}: {manifests}"
                )
            for manifest, outputs in missing_outputs.items():
                logger.error(f"Manifest {manifest} lists {len(outputs)} missing output files: {outputs}")
            return False

        total_expected = sum(len(v) for v in expected_by_biobank.values())
        logger.success(f"Verified {total_expected} STEP2 block output lists and listed files")
        return True

    # ------------------------------------------------------------------
    # Meta-analysis job
    # ------------------------------------------------------------------

    async def _build_meta_job(self, meta_cfg: MetaAnalysisConfig) -> list[Job]:
        """Build the meta-analysis jobs: one applet job per chromosome for `local`, else one WDL job.

        Raises:
            ValueError: If `meta_workflow` is unset.
        """
        if not meta_cfg.meta_workflow:
            raise ValueError("meta_analysis_config.meta_workflow must be set for meta-analysis")

        if meta_cfg.location == RunEnv.LOCAL:
            jobs: list[Job] = []
            for chrom in self.config.chromosomes_to_run:
                job = Job(
                    biobank=BioBank.UKB,
                    description=f"{self.config.run_prefix}_meta_analysis_chr{chrom}",
                    recipe={"op": "meta", "chr": chrom},
                )
                job._parent = weakref.ref(self)
                await job.generate_input()
                jobs.append(job)
            return jobs
        else:
            job = Job(
                biobank=BioBank.AOU,
                description=f"{self.config.run_prefix}_meta_analysis_discovery",
                recipe={"op": "meta"},
            )
            job._parent = weakref.ref(self)
            await job.generate_input()
            return [job]

    # ------------------------------------------------------------------
    # Step-0 subset
    # ------------------------------------------------------------------

    def _meta_biobank(self) -> BioBank:
        """Return the biobank whose platform runs the meta and step-0 subset jobs."""
        meta_cfg = self.config.meta_analysis_config
        return BioBank.UKB if meta_cfg is not None and meta_cfg.location == RunEnv.LOCAL else BioBank.AOU

    def _meta_output_folder(self) -> str:
        """Return the meta job's output folder: a GCS path without `gs://`, or a DNAnexus path."""
        if self._meta_biobank() is BioBank.UKB:
            return f"{self.config.cohorts[BioBank.UKB].output_location}/meta_analysis/all_biobank_meta"
        return f"{self._meta_base_folder_gcs()}/all_biobank_meta"

    def _step0_folder(self, cond_type: CondType) -> str:
        """Return the folder of a type's per-chromosome step-0 subsets at the configured threshold."""
        return f"{self._meta_output_folder()}/{self.config.step0_subset_dirname(cond_type)}"

    def _step0_cond_types(self) -> list[CondType]:
        """Return the conditioning types with step-0 subsets, from the stages the meta mode runs."""
        stages = self.meta_mode.stages(self.config.has_cond_file)
        return [
            cond_type
            for cond_type, needed in (
                (CondType.SINGLEVAR, stages.uncond_singles),
                (CondType.AGGREGATE, stages.aggregates),
            )
            if needed
        ]

    async def _missing_step0_chromosomes(self, folder: str) -> list[str]:
        """Return the chromosomes without a completed subset in `folder`, judged by the task logs."""
        prefix = re.escape(self.config.run_prefix)
        try:
            done = set(
                await self.platforms[self._meta_biobank()].list_files(
                    folder=folder, pattern=rf"^{prefix}_step0_subset_chr.+\.log$"
                )
            )
        except Exception as exc:  # noqa: BLE001 — a folder not created yet is an error on DNAnexus
            logger.debug(f"No completed step-0 subsets listed in {folder}: {exc}")
            done = set()
        return [
            c
            for c in self.config.chromosomes_to_run
            if f"{self.config.run_prefix}_step0_subset_chr{c}.log" not in done
        ]

    async def _build_step0_jobs(self) -> list[Job]:
        """Build subset jobs for the chromosomes not yet subset at the configured thresholds."""
        meta_bb = self._meta_biobank()
        jobs: list[Job] = []
        for cond_type in self._step0_cond_types():
            folder = self._step0_folder(cond_type)
            missing = await self._missing_step0_chromosomes(folder)
            if not missing:
                logger.info(f"{cond_type.value} step-0 subsets already complete in {folder}")
                continue
            logger.info(f"Subsetting {cond_type.value} results for {len(missing)} chromosome(s) into {folder}")
            # A WDL scatters over chromosomes itself; an applet job covers one
            chunks = [missing] if meta_bb is BioBank.AOU else [[c] for c in missing]
            for chroms in chunks:
                suffix = f"_chr{chroms[0]}" if meta_bb is BioBank.UKB else ""
                job = Job(
                    biobank=meta_bb,
                    description=f"{self.config.run_prefix}_step0_subset_{cond_type.value}{suffix}",
                    recipe={
                        "op": "step0_subset",
                        "cond_type": cond_type.value,
                        "chromosomes": chroms,
                        "output_folder": folder,
                    },
                )
                job._parent = weakref.ref(self)
                await job.generate_input()
                jobs.append(job)
        return jobs

    async def _run_step0_subset(
        self, reset_jobs: JobState | None, reset_jobs_list: list[str], attempt_job_reset: bool
    ) -> None:
        """Build any missing step-0 subsets remotely, then merge them into the run directory."""
        if self.config.meta_analysis_config is None:
            logger.error("meta_analysis_config must be set in the config to run the step0_subset phase")
            self.had_error = True
            return
        folders = {t: self._step0_folder(t) for t in self._step0_cond_types()}
        # Earlier jobs are resumed only while unfinished and for the current thresholds;
        # otherwise they are rebuilt for the chromosomes still missing
        if self.jobs and (
            all(j.state == JobState.COMPLETED for j in self.jobs)
            or any(j.recipe.get("output_folder") not in folders.values() for j in self.jobs)
        ):
            self.jobs.clear()
        await self._run_jobs(self._build_jobs_for_phase, reset_jobs, reset_jobs_list, attempt_job_reset)
        if self.had_error:
            return

        incomplete = {t: m for t, f in folders.items() if (m := await self._missing_step0_chromosomes(f))}
        if incomplete:
            for cond_type, missing in incomplete.items():
                logger.error(
                    f"No completed {cond_type.value} step-0 subset for chromosome(s) {', '.join(missing)} "
                    f"in {folders[cond_type]}"
                )
            self.had_error = True
            return

        self._warn_started_conditioning(list(folders))
        await self._fetch_and_merge_step0(folders)
        self.persist_state(reason="step0_merge_complete", max_delay=0)

    def _warn_started_conditioning(self, cond_types: list[CondType]) -> None:
        """Warn for phenotypes whose conditioning has started, as it keeps its existing step-0 subset."""
        from fedcat.core.conditioning import ConditioningStep

        if not self.config.cond_suffix:
            return
        for cond_type in cond_types:
            for pheno in self.config.pheno_names:
                step = ConditioningStep.from_discovery_config(self.config, pheno, cond_type)
                if (step.run_path / f"{cond_type.value}_cond_step0").exists():
                    logger.warning(
                        f"{cond_type.value} conditioning for {pheno!r} has already started in {step.run_path}; "
                        "it keeps using the step-0 subset it started from. Move that directory aside to "
                        "restart conditioning from the new subset."
                    )

    async def _fetch_and_merge_step0(self, folders: dict[CondType, str]) -> None:
        """Download each phenotype's per-chromosome subsets and merge them into the run directory.

        Args:
            folders: Subset folder per conditioning type.
        """
        platform = self.platforms[self._meta_biobank()]
        merge_root = self.run_path / ".meta_step0_merged"
        for pheno in self.config.pheno_names:
            with log_context(PHENO_VAR, pheno):
                pheno_dl = merge_root / f"{pheno}_parts"
                if pheno_dl.exists():
                    shutil.rmtree(pheno_dl)
                pheno_dl.mkdir(parents=True)
                for cond_type, folder in folders.items():
                    # Removed even when no subset comes back, so a stale merge is never reused
                    (self.run_path / f"{pheno}_{cond_type.value}_cond_step0_subset.tar.gz").unlink(
                        missing_ok=True
                    )
                    await platform.download(
                        folder=folder,
                        pattern=rf"^{re.escape(pheno)}_{cond_type.value}_cond_step0_subset_chr[^_]+\.tar\.gz$",
                        dest=pheno_dl,
                        max_workers=4,
                    )
                await asyncio.to_thread(_merge_step0_subsets, pheno, pheno_dl, self.run_path)
        shutil.rmtree(merge_root, ignore_errors=True)

    # ------------------------------------------------------------------
    # Conditioning
    # ------------------------------------------------------------------

    @staticmethod
    def _retry_failed_conditioning(step) -> None:
        """Apply `--force-job-retry` to a conditioning step.

        Failed jobs are resubmitted; otherwise a failed step reruns the stage it failed in,
        and a step without failures resumes as it was.
        """
        from fedcat.core.structures import StepState

        failed = [j for j in step.jobs if j.state == JobState.ERROR]
        if failed:
            step.set_state(StepState.ERROR)
            step.reset_error_state(target_state=StepState.RUN_JOBS, target_job_state=JobState.EMPTY)
            logger.info(f"Retrying {len(failed)} failed job(s)")
        elif step.state == StepState.ERROR:
            step.reset_error_state()
            logger.info(f"Retrying from {step.state.name}")
        else:
            logger.info("No failed jobs to retry; resuming")

    @staticmethod
    def _autodetect_cond_step(
        run_path: Path,
        cond_type: CondType,
    ) -> int:
        """Return the highest `N` of the `{cond_type}_cond_step{N}` directories in `run_path`, or 0."""
        if not run_path.is_dir():
            return 0
        pattern = rf"^{re.escape(cond_type.value)}_cond_step(\d+)$"
        steps = [
            int(m.group(1))
            for item in run_path.iterdir()
            if item.is_dir() and (m := re.match(pattern, item.name))
        ]
        if steps:
            detected = max(steps)
            logger.info(f"Auto-detected most recent conditioning step {detected} in {run_path}")
            return detected
        return 0

    async def _run_conditioning(self, cond_type: CondType) -> CondOutcome:
        """Run forward-selection conditioning for every phenotype not yet completed, concurrently.

        Each phenotype resumes from its saved step state (the step given by `--cond-step`,
        or the latest step directory) or starts at step 0 from its step-0 subset. Steps
        run until no signals remain, `cond_max_steps` is reached, or `--cond-one-step`
        stops after one. A saved config that differs from the current one is an error
        unless `--overwrite-config` is given.

        Args:
            cond_type: Conditioning type.

        Returns:
            FAILED if any phenotype failed, else INCOMPLETE if any stopped early, else COMPLETED.
        """
        # Deferred: conditioning imports this module
        from fedcat.core.conditioning import ConditioningStep
        from fedcat.core.structures import StepState

        base_path = Path(self.config.local_run_path)
        base_path.mkdir(parents=True, exist_ok=True)

        phenos_to_run = [
            ph
            for ph in self.config.pheno_names
            if self.cond_status.get(cond_type, {}).get(ph) is not CondOutcome.COMPLETED
        ]

        max_steps = self.config.cond_max_steps
        one_step = self.cond_one_step
        overwrite_config = self._overwrite_config
        conv = make_conv()

        # A step-0 tarball given on the command line holds one phenotype's subset
        override_tarball = self.cond_singles_tarball if cond_type == CondType.SINGLEVAR else self.cond_agg_tarball
        if override_tarball is not None and len(phenos_to_run) > 1:
            logger.error(
                f"Direct step-0 tarball input for {cond_type.value} conditioning "
                f"({override_tarball}) is only supported for a single-phenotype run, "
                f"but {len(phenos_to_run)} phenotypes are pending: {', '.join(phenos_to_run)}."
            )
            return CondOutcome.FAILED

        async def _run_one_pheno(pheno: str) -> CondOutcome:
            """Resume or start one phenotype's conditioning, and run it to an outcome."""
            PHENO_VAR.set(pheno)
            # 1. Step to start from; the probe step locates the phenotype's directory
            probe = ConditioningStep.from_discovery_config(
                config=self.config,
                pheno_name=pheno,
                cond_type=cond_type,
            )
            pheno_run_path = probe.run_path

            if self.cond_step is not None:
                effective_step = self.cond_step
            else:
                effective_step = self._autodetect_cond_step(pheno_run_path, cond_type)

            # 2. Load the saved step, or create step 0
            probe.step_number = effective_step
            state_file = probe.state_path

            if state_file.exists():
                logger.info(f"Loading saved conditioning state from {state_file}")
                raw_step = json.loads(state_file.read_text())
                for key in pin_legacy_config(raw_step.get("config", {})):
                    logger.warning(
                        f"{state_file} predates {key!r}; keeping its earlier behaviour ({key}=false). "
                        f"Set it to false in the config to resume, or restart conditioning from the "
                        f"meta-analysis to use true."
                    )
                step = conv.structure(raw_step, ConditioningStep)
                if step.step_number != effective_step or step.cond_type != cond_type:
                    logger.error(
                        f"{state_file} holds {step.cond_type.value} step {step.step_number} state, not "
                        f"{cond_type.value} step {effective_step}; the file does not match its name, so fix "
                        f"or remove it rather than resume from it."
                    )
                    return CondOutcome.FAILED

                diffs = config_differences(self.config, step.config)
                adopt_run_controls(step.config, self.config)
                if diffs and overwrite_config:
                    step.adopt_config(self.config)
                    logger.info(
                        f"Updated the saved conditioning config (--overwrite-config):\n"
                        f"{format_config_differences(diffs)}"
                    )
                elif diffs:
                    logger.error(
                        f"The conditioning state was saved with a different config:\n"
                        f"{format_config_differences(diffs)}\nUse --overwrite-config to replace it."
                    )
                    return CondOutcome.FAILED
            else:
                if effective_step > 0:
                    logger.error(
                        f"No saved conditioning state for step {effective_step} "
                        f"({state_file}), but step directories suggest progress beyond step 0. "
                        f"To fix this, delete or move the problematic step directories."
                    )
                    return CondOutcome.FAILED
                logger.info(f"No saved conditioning state; starting fresh at step {effective_step}")
                step = ConditioningStep.from_discovery_config(
                    config=self.config,
                    pheno_name=pheno,
                    cond_type=cond_type,
                    tarball_path=(
                        override_tarball
                        or self.run_path.resolve() / f"{pheno}_{cond_type.value}_cond_step0_subset.tar.gz"
                    ),
                    step_number=effective_step,
                )

                # A later step directory means the wrong step was chosen to resume from
                next_step_dir = pheno_run_path / f"{cond_type.value}_cond_step{effective_step + 1}"
                if next_step_dir.exists():
                    logger.error(
                        f"Step {effective_step + 1} directory already exists "
                        f"({next_step_dir}).  Are you sure you chose the "
                        f"correct step to resume?  If so, delete or move the "
                        f"directory."
                    )
                    return CondOutcome.FAILED

            COND_STEP_VAR.set(step)

            # 3. Apply the reset flags
            if self.cond_retry_failed:
                self._retry_failed_conditioning(step)
            elif self.cond_reset_step is not None or self.cond_reset_jobs is not None:
                target_state = StepState(self.cond_reset_step) if self.cond_reset_step is not None else None
                target_job_state = JobState(self.cond_reset_jobs) if self.cond_reset_jobs is not None else None
                if step.state == StepState.ERROR:
                    step.reset_error_state(
                        target_state=target_state,
                        target_job_state=target_job_state,
                    )
                    logger.info(
                        "Reset conditioning step from ERROR"
                        + (f" to {target_state.name}" if target_state else "")
                        + (f" (jobs → {target_job_state.name})" if target_job_state else "")
                    )
                elif self.cond_force_reset:
                    step.set_state(StepState.ERROR)
                    step.reset_error_state(
                        target_state=target_state,
                        target_job_state=target_job_state,
                    )
                    logger.warning(
                        "Force-reset conditioning step"
                        + (f" to {target_state.name}" if target_state else "")
                        + (f" (jobs → {target_job_state.name})" if target_job_state else "")
                        + " (--force-reset)"
                    )
                else:
                    logger.error(
                        f"--cond-reset-step/--reset-jobs specified but the conditioning step "
                        f"is not in ERROR state ({step.state.name}).  "
                        f"Use --force-reset to override."
                    )
                    return CondOutcome.FAILED

            # 4. Run steps until done
            result = CondOutcome.FAILED
            try:
                while True:
                    success = await step.run_to_completion()
                    if not success:
                        logger.error("Conditioning step failed")
                        result = CondOutcome.FAILED
                        break
                    if not step.required:
                        if await step.write_final_results():
                            result = CondOutcome.COMPLETED
                        else:
                            logger.error("Could not write the final conditioning results; resume to retry")
                        break
                    if one_step:
                        logger.info("--cond-one-step: stopping after this step")
                        result = CondOutcome.INCOMPLETE
                        break
                    if 0 < max_steps <= step.step_number:
                        logger.warning(
                            f"Reached cond_max_steps={max_steps} with signals left; writing partial results "
                            f"(raise cond_max_steps and resume with --start-phase {self.current_phase.value} "
                            "to carry on)"
                        )
                        if await step.write_final_results(partial=True):
                            result = CondOutcome.INCOMPLETE
                        else:
                            logger.error("Could not write the partial conditioning results; resume to retry")
                        break
                    # The state file is per step, so save this step's before moving on
                    await step.flush_state()
                    try:
                        step.increment()
                    except (FileExistsError, RuntimeError) as exc:
                        logger.error(f"increment() failed: {exc}")
                        result = CondOutcome.FAILED
                        break
            except asyncio.CancelledError:
                logger.info("Saving conditioning state before pausing…")
                raise
            except Exception:
                step.set_state(StepState.ERROR)
                logger.exception("Unhandled exception in conditioning loop")
                result = CondOutcome.FAILED
            finally:
                step.persist_state(reason="loop_exit", max_delay=0)
                await step.flush_state()
            return result

        # Phenotypes start 10 s apart, to spread their submissions
        async def _staggered(pheno: str, delay: float) -> CondOutcome:
            """Run one phenotype after `delay` seconds."""
            if delay > 0:
                await asyncio.sleep(delay)
            return await _run_one_pheno(pheno)

        if not phenos_to_run:
            logger.info("All phenotypes already complete; nothing to do")
        else:
            logger.info(
                f"Starting conditioning for {len(phenos_to_run)} phenotype(s): {', '.join(phenos_to_run)}"
            )

        tasks: list[asyncio.Task] = [
            asyncio.create_task(_staggered(pheno, i * 10.0)) for i, pheno in enumerate(phenos_to_run)
        ]

        logger.trace(
            f"Launching conditioning tasks for phenotypes: {phenos_to_run} with delays: {[i * 10.0 for i in range(len(phenos_to_run))]}"
        )
        try:
            results: list[CondOutcome | BaseException] = await asyncio.gather(*tasks, return_exceptions=True)
        except asyncio.CancelledError:
            # Let every phenotype save its state, then pass the pause on: a paused
            # phenotype has not failed, so it must not be recorded as one.
            await asyncio.shield(asyncio.gather(*tasks, return_exceptions=True))
            raise

        logger.trace(f"conditioning task results: {results}")

        # Only a completed phenotype is skipped next time; incomplete or failed ones resume
        status = self.cond_status.setdefault(cond_type, {})
        for pheno, res in zip(phenos_to_run, results):
            if isinstance(res, BaseException):
                logger.error(f"Conditioning task for {pheno!r} raised: {res!r}")
                res = CondOutcome.FAILED
            status[pheno] = res

        logger.trace(f"cond_status so far: {self.cond_status}")
        self.persist_state(reason=f"conditioning_{cond_type.value}_status_updated", max_delay=0)
        for outcome in (CondOutcome.FAILED, CondOutcome.INCOMPLETE):
            if phenos := [ph for ph in phenos_to_run if status[ph] is outcome]:
                level = "ERROR" if outcome is CondOutcome.FAILED else "WARNING"
                logger.log(level, f"{cond_type.value} conditioning {outcome.value} for: {', '.join(phenos)}")
                return outcome
        return CondOutcome.COMPLETED
