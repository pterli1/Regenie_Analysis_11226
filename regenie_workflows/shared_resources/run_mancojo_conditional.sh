#!/bin/bash
# Derive conditional summary statistics with manc-cojo, as an alternative to
# re-running regenie with --condition-list.
#
# Why this exists: regenie re-fits a joint logistic model at every conditioning
# step, and with a large block of ultra-rare conditioning covariates on an
# imbalanced binary trait that model separates, giving non-convergence and
# inflated effect estimates. COJO never fits a joint model — it projects marginal
# estimates through an LD matrix — so separation is structurally impossible. The
# trade is that the projection is a linear-model approximation applied to
# log-odds effects, which degrades as odds ratios and imbalance grow.
#
# Because COJO works from marginal statistics, the .ma passed here must come from
# an UNCONDITIONAL regenie run over the same samples and variants as the pgen:
# COJO assumes the LD reference is the GWAS sample itself.
set -euxo pipefail

pgen=${1}
ma_file=${2}
conditional_list=${3}
chr=${4}
output_prefix=${5}
threads=${6}
collinear=${7:-0.9}
# --effect-mode: GCTA, imputeNA or removeNA (see cond_mancojo_effect_mode in fedcat)
effect_mode=${8:-imputeNA}
# Estimate genotypic variance from the reference rather than assuming 2pq. For
# single variants under HWE the two agree closely; for a collapsed burden mask
# 2pq is a poor estimate, so this matters once masks enter the .ma.
var_from_ld=${9:-1}

# The WDL/applet stages dependencies into the task working directory rather than
# onto PATH, so prefer ./manc_cojo and fall back to an installed one.
if [[ -x ./manc_cojo ]]; then
  manc_cojo_bin=./manc_cojo
elif command -v manc_cojo >/dev/null 2>&1; then
  manc_cojo_bin=manc_cojo
else
  echo "FATAL: manc_cojo not found in the task directory or on PATH." >&2
  echo "       Add the binary at regenie_workflows/shared_resources/manc_cojo so" >&2
  echo "       'fedcat setup-workflows' stages it with the other dependencies." >&2
  exit 1
fi

com="${manc_cojo_bin} \
  --pfile ${pgen} \
  --cojo-file ${ma_file} \
  --cojo-cond ${conditional_list} \
  --cojo-collinear ${collinear} \
  --effect-mode ${effect_mode} \
  --chr ${chr} \
  --thread-num ${threads} \
  --out ${output_prefix}"

if [ "${var_from_ld}" == "1" ]; then
  com="${com} \
  --var-from-ld"
fi

$com

# An absent .cma.cojo here should mean every candidate was conditioned away or screened out
if [[ ! -s "${output_prefix}.cma.cojo" ]]; then
  echo "[manc-cojo] no candidates left for chr${chr}" >&2
fi
