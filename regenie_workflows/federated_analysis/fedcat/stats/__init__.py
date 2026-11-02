"""Statistics: meta-analysis, ACAT, collinearity and conditioning inputs and summaries.

- [`meta_analysis`][fedcat.stats.meta_analysis]: IVW fixed-effects meta-analysis of
  REGENIE results (`fedcat meta`).
- [`acat`][fedcat.stats.acat]: ACAT tests of aggregates and annotation conversion
  (`fedcat acat`, `fedcat convert-annots`).
- [`subset_results`][fedcat.stats.subset_results]: step-0 subsets of significant
  results for conditioning (`fedcat subset-results`).
- [`collinearity`][fedcat.stats.collinearity]: collinearity of candidate variants with
  the selected set (`fedcat collinearity`).
- [`summarise_indep`][fedcat.stats.summarise_indep]: conditioning summaries of finished
  runs (`fedcat summarise-indep`).
"""

from .acat import compute_acat
from .collinearity import run_collinearity_check
from .meta_analysis import (
    log10p_to_p,
    meta_analyse_regenie,
    meta_analyse_regenie_async,
    z_to_log10p,
)

__all__ = [
    "compute_acat",
    "log10p_to_p",
    "meta_analyse_regenie",
    "meta_analyse_regenie_async",
    "run_collinearity_check",
    "z_to_log10p",
]
