"""`fedcat summarise-indep`: rewrite the conditioning results of a finished run.

Uses only what every fedcat version leaves on disk: each step's `summary_output`, and
the `passing_ids` key of a step state file. The state file is read as plain JSON
rather than as a `ConditioningStep`, so state files of older versions also work.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import polars as pl
from loguru import logger

from fedcat.core.conditioning import generate_final_output
from fedcat.core.structures import CondType

# Step number of a step directory
_STEP_RE = re.compile(r"_cond_step(\d+)$")


def _detect_end_step(run_dir: Path, cond_type: CondType) -> int:
    """Return the last step of a run with a `summary_output`.

    A step that ended before its meta-analysis has none, so results end at the step
    before, as when conditioning writes them.

    Raises:
        FileNotFoundError: If no step has a `summary_output`.
    """
    steps: list[int] = []
    for d in run_dir.glob(f"{cond_type.value}_cond_step*"):
        if not d.is_dir():
            continue
        m = _STEP_RE.search(d.name)
        if m and (d / "summary_output").is_file():
            steps.append(int(m.group(1)))
    if not steps:
        raise FileNotFoundError(
            f"No '{cond_type.value}_cond_step*' directories with a summary_output found under {run_dir}"
        )
    return max(steps)


def _state_json_templates(cond_type: CondType) -> list[tuple[str, re.Pattern]]:
    """Return the `(glob, step_number_regex)` of the state file names, current first.

    Older versions wrote `step{N}_config.json` rather than
    `conditioning_state_{cond}_step{N}.json`.
    """
    return [
        (
            f"conditioning_state_{cond_type.value}_step*.json",
            re.compile(rf"conditioning_state_{re.escape(cond_type.value)}_step(\d+)\.json$"),
        ),
        (
            "step*_config.json",
            re.compile(r"step(\d+)_config\.json$"),
        ),
    ]


def _locate_passing_ids_json(run_dir: Path, cond_type: CondType, end_step: int) -> Path:
    """Return the state file of step `end_step`, else of the highest step found.

    Raises:
        FileNotFoundError: If there is no state file.
    """
    candidates: list[tuple[int, Path]] = []
    for glob_pat, step_re in _state_json_templates(cond_type):
        for p in run_dir.glob(f"**/{glob_pat}"):
            m = step_re.search(p.name)
            if m:
                candidates.append((int(m.group(1)), p))

    if not candidates:
        raise FileNotFoundError(
            f"Could not locate a 'conditioning_state_{cond_type.value}_step*.json' or "
            f"'step*_config.json' under {run_dir}; pass --passing-ids-json explicitly."
        )

    for n, p in candidates:
        if n == end_step:
            return p
    return max(candidates, key=lambda np: np[0])[1]


def _detect_stat_col(run_dir: Path, cond_type: CondType, end_step: int) -> str:
    """Return `"LOG10P"`, or `"P"` for older runs, from the columns of the last `summary_output`.

    Raises:
        ValueError: If it has neither column.
        FileNotFoundError: If no step has a `summary_output`.
    """
    for i in range(end_step, -1, -1):
        path = run_dir / f"{cond_type.value}_cond_step{i}" / "summary_output"
        if not path.is_file():
            continue
        cols = pl.read_csv(path, separator="\t", n_rows=1, ignore_errors=True).columns
        if "LOG10P_OUT" in cols or "LOG10P_IN" in cols:
            return "LOG10P"
        if "P_OUT" in cols or "P_IN" in cols:
            return "P"
        raise ValueError(f"No LOG10P/P test-statistic columns found in {path}; columns: {cols}")
    raise FileNotFoundError(f"No 'summary_output' file found under {run_dir} for steps 0..{end_step}")


def _read_passing_ids(json_path: Path) -> dict[str, list[str]]:
    """Read `passing_ids` from a state file, as string chromosomes and IDs.

    Raises:
        ValueError: If the file has no `passing_ids` mapping of lists.
    """
    with open(json_path) as fh:
        data = json.load(fh)
    if not isinstance(data, dict) or "passing_ids" not in data:
        raise ValueError(f"No top-level 'passing_ids' key in {json_path}")
    passing_ids = data["passing_ids"]
    if not isinstance(passing_ids, dict) or not all(isinstance(v, list) for v in passing_ids.values()):
        raise ValueError(f"'passing_ids' in {json_path} is not a mapping of chrom -> list of IDs")
    return {str(chrom): [str(i) for i in ids] for chrom, ids in passing_ids.items()}


def summarise_independent_results(
    *,
    run_dir: Path,
    cond_type: CondType,
    pheno_name: str,
    passing_ids_json: Path | None = None,
    end_step: int | None = None,
    output_dir: Path | None = None,
    stat_col: str | None = None,
) -> Path | None:
    """Write a finished conditioning run's results; see [`generate_final_output`][fedcat.core.conditioning.generate_final_output].

    Args:
        run_dir: Conditioning run directory holding the step directories.
        cond_type: Conditioning type.
        pheno_name: Phenotype name, used in the output filenames.
        passing_ids_json: State file holding `passing_ids`; found in `run_dir` when None.
        end_step: Last step to include; the last with a `summary_output` when None.
        output_dir: Results directory; `<type>_cond_results` in `run_dir` when None.
        stat_col: `"LOG10P"`, or `"P"` for older runs; detected when None.

    Returns:
        The results directory, or None if nothing was written (e.g. nothing was
        selected at step 0, or the results directory exists).

    Raises:
        FileNotFoundError: If `run_dir`, a summary or a state file is missing.
    """
    if not run_dir.is_dir():
        raise FileNotFoundError(f"Run directory not found: {run_dir}")

    if end_step is None:
        end_step = _detect_end_step(run_dir, cond_type)
        logger.info(f"Auto-detected end step: {end_step}")

    if passing_ids_json is None:
        passing_ids_json = _locate_passing_ids_json(run_dir, cond_type, end_step)
        logger.info(f"Using passing_ids JSON: {passing_ids_json}")

    if stat_col is None:
        stat_col = _detect_stat_col(run_dir, cond_type, end_step)
        logger.info(f"Auto-detected test-statistic column: {stat_col}")

    passing_ids = _read_passing_ids(passing_ids_json)
    n_passing = sum(len(v) for v in passing_ids.values())
    logger.info(f"Loaded {n_passing} passing IDs across {len(passing_ids)} chromosomes")

    results_dir = generate_final_output(
        run_path=run_dir,
        cond_type=cond_type,
        passing_ids=passing_ids,
        pheno_name=pheno_name,
        required=False,
        end_step=end_step,
        output_dir=output_dir,
        stat_col=stat_col,
    )
    if results_dir is None:
        logger.warning("No summary outputs were produced")
    else:
        logger.success(f"Summaries written to: {results_dir.resolve()}")
    return results_dir
