#!/bin/bash
set -euxo pipefail

run_prefix=${1}
catvars=${2}
output_directory=${3}
output_prefix=${4}
chr=${5}
ncon=${6}
masks=${7}
annotations=${8}
setlist=${9}
mafs=${10}
invnorm=${11}
runtype=${12}
pheno_type=${13}
exclvars=${14}
phenos_to_run=${15}
threads=${16}
# Write the burden masks to PLINK bed as well as testing them. The manc-cojo
# engine needs the mask pseudo-genotypes to compute LD between masks and the
# single variants being conditioned on; regenie builds them on the fly and this
# is the only way to get them out. Requires the default --build-mask 'max'
# ('sum' is documented as incompatible with --write-mask).
write_mask=${17:-0}

# create output filename set
if [ ${runtype} == "mac10" ];then
  output_name="${output_directory}/${output_prefix}_mac10"
elif [ ${runtype} == "full" ];then
  output_name="${output_directory}/${output_prefix}_burden"
fi

# build main command
com="regenie  \
  --step 2 \
  --phenoFile phenotype_tsv \
  --covarFile covariate_tsv \
  --bsize 400 \
  --minMAC=5 \
  --max-condition-vars 20000 \
  --par-region b38 \
  --anno-file ${annotations} \
  --set-list ${setlist} \
  --mask-def ${masks} \
  --threads ${threads} \
  --pred ${run_prefix}_Step1_pred.list \
  --aaf-bins ${mafs} \
  --out ${output_name} \
  --phenoColList ${phenos_to_run} \
  --gz"

# add other args as required
# which subset pgen to use
if [ ${runtype} == "mac10" ]; then
  com=${com}" \
	--pgen chr${chr}_for_regenie_mac10"
elif [ ${runtype} == "full" ]; then
	com=${com}" \
	--pgen chr${chr}_for_regenie"
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
  --condition-list conditional_list_local \
  --condition-file pgen,chr${chr}_for_conditioning"
fi

# if mask genotypes are needed downstream
if [ "${write_mask}" == "1" ]; then
  com=${com}" \
  --write-mask"
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

# if covariates to exclude
if [ "${exclvars}" != "None" ]; then
  com=${com}" \
  --covarExcludeList ${exclvars}"
fi

$com
