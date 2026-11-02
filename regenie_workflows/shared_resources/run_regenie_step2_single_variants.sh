#!/bin/bash
set -euxo pipefail

run_prefix=${1}
catvars=${2}
output=${3}
chr=${4}
invnorm=${5}
ncon=${6}
output_prefix=${7}
pgen=${8}
cond_pgen=${9}
pheno_type=${10}
exclvars=${11}
phenos_to_run=${12}
threads=${13}
interaction=${14:-}
# REGENIE's own default is 5. The manc-cojo engine needs marginal statistics for
# the ultra-rare variants it conditions on, which are almost all below that, so
# it passes 1; every other caller leaves this unset and the flag is not emitted.
minmac=${15:-}

# build main command
com="regenie  \
  --step 2 \
  --pgen ${pgen} \
  --out ${output}/${output_prefix} \
  --phenoFile phenotype_tsv \
  --covarFile covariate_tsv \
  --bsize 400 \
  --max-condition-vars 20000 \
  --par-region b38 \
  --threads ${threads} \
  --pred ${run_prefix}_Step1_pred.list \
  --phenoColList ${phenos_to_run} \
  --gz"

# add other args as required
# if a non-default minor allele count floor was asked for
if [ -n "${minmac}" ]; then
  com=${com}" \
  --minMAC=${minmac}"
fi

# if categorical covariates
if [[ ${catvars} != "NONE" && ${catvars} != "" ]]; then
	com=${com}" \
	--catCovarList ${catvars} \
	--maxCatLevels 30"
fi

# if transform has been asked for
if [ "${invnorm}" == "RINT" ];then
  com=${com}" \
  --apply-rint"
elif [ "${invnorm}" == "MCC" ]; then
  com=${com}" \
  --mcc \
  --mcc-thr 0.01"
fi

# if there are variants to condition on
if [ $ncon -ne 0 ]; then
  com=${com}" \
  --condition-list conditional_list_local"
  # if should use cond pgen
  if [ $cond_pgen != "None" ]; then
    com=${com}" \
    --condition-file pgen,${cond_pgen}"
  fi
fi

# if binary pheno
if [ "${pheno_type}" == "Binary" ];then
  com=${com}" \
  --bt \
  --firth \
  --firth-se \
  --approx \
  --pThresh 0.01"
fi

# if testing an interaction (variable must be a column of the covariate file)
if [ "${interaction}" != "" ]; then
  com=${com}" \
  --interaction ${interaction}"
fi

# if covariates to exclude
if [ "${exclvars}" != "None" ]; then
  com=${com}" \
  --covarExcludeList ${exclvars}"
fi

# run regenie
$com
