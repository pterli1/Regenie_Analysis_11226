# Shared Resources

This folder contains binaries, scripts and other shared resources referenced by multiple workflows in this repository.

## Contents
- `plink2` — plink2 binary (ignored by git, must be added manually)
- `regenie_v<version>.gz_x86_64_ubuntu20_mkl` — REGENIE binaries for the UKB applets, one per supported version (e.g. `regenie_v4.1.3...`, `regenie_v4.1...`; ignored by git, must be added manually). Jobs pick one with the config's `regenie_version`; the default `4.1.3` must be present
- `manc_cojo` — Manc-COJO binary, Linux x86_64 (ignored by git, must be added manually). Only needed for
  `cond_engine = "mancojo"`; build it from the Exeter fork, which adds pgen input and the chromosome
  parsing REGENIE uses. It is the one optional dependency: `fedcat setup-workflows` stages it when
  present and warns when it is not, and its absence is only a pre-run error once `cond_engine` is set
  to `mancojo`.
- `run_regenie_step2_single_variants.sh` - executes REGENIE step 2 binary; translates positional inputs into REGENIE args with a number of defaults args suitable for single variant analysis
- `run_regenie_step2_aggregates.sh` - executes REGENIE step 2 binary; translates positional inputs into REGENIE args with a number of defaults args suitable for aggregate analysis
- `run_mancojo_conditional.sh` - derives conditional summary statistics with manc-cojo instead of re-running
  REGENIE, for binary traits where the joint logistic model separates on a large block of rare conditioning
  covariates. Requires an unconditional `.ma` built by `build_mancojo_ma` (in `shared_functions.sh`).

## Notes
- These files are intended to be used by the REGENIE-based WDL workflows and DNAnexus applets
