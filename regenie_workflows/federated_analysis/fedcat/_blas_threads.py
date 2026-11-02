"""Set BLAS threads from `fedcat collinearity --threads` on import, before numpy loads.

BLAS reads its thread count when numpy is first imported, which the CLI does before
parsing arguments, so `__main__` imports this module first.
"""

import os
import sys

BLAS_THREAD_VARS = (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)


def preset_blas_threads(argv: list[str]) -> None:
    """Set the BLAS thread environment variables from `--threads` when `argv` runs `collinearity`.

    Args:
        argv: Command-line arguments, as `sys.argv`.
    """
    if argv[1:2] != ["collinearity"]:
        return
    threads = None
    for i, arg in enumerate(argv):
        if arg in ("-t", "--threads") and i + 1 < len(argv):
            threads = argv[i + 1]
        elif arg.startswith("--threads="):
            threads = arg.split("=", 1)[1]
    if threads:
        for var in BLAS_THREAD_VARS:
            os.environ[var] = threads


preset_blas_threads(sys.argv)
