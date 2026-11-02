#!/bin/bash
set -euxo pipefail

phenotype_tsv=${1}
covariate_tsv=${2}
exclcovars=${3}
conditional=${4}
type=${5}
invnorm=${6}
catvars=${7}
chr=${8}
vars_to_run=${9}
condition_vars=${10}
stratum=${11}
bucket=${12}
run_prefix=${13}
run_suffix=${14}
pheno_to_run=${15}
check_coll=${16}
step1_directory=${17}
step1_input_template=${18}
results_directory=${19}
pfile_template=${20}
# Prefix under which step1 predictions were written; defaults to run_prefix so
# older callers that don't pass it keep working.
step1_run_prefix=${21:-$run_prefix}
# Covariate stratification: covariate column and the value to keep. Empty (the default
# for callers that don't pass them) means run on everyone.
stratify_covar=${22:-}
stratify_value=${23:-}
# Engine used to derive conditional statistics: "regenie" re-fits an adjusted model
# per step; "mancojo" projects the step0 marginal statistics through an LD matrix.
cond_engine=${24:-regenie}
cond_rsq_threshold=${25:-0.9}
# manc-cojo --effect-mode (GCTA, imputeNA or removeNA)
cond_mancojo_effect_mode=${26:-imputeNA}
# Conditioning step number, from fedcat
cond_step=${27}

source shared_functions.sh

main() {
  version="1.7.1"
  # resolve step1 template into step1_input_path/_directory/_name
  split_template "${step1_input_template}" stratum "${stratum}" step1_input
  # process selected vars to condition on
  if [ -z "$condition_vars" ]; then
    touch additional_vars_conditioning
  else
    echo "${condition_vars}" | sed 's/\,/\n/g' > additional_vars_conditioning
  fi
  step=${cond_step}
  job_prefix=${pheno_to_run}_${stratum}_chr${chr}_singlevar_conditioning_step${step}
  log_file="${job_prefix}.workflow.log"
  threads=$(nproc)
  # resolve pfile template into pfile_in_path/_directory/_name (pfile_name below is a derived local name)
  split_template "${pfile_template}" chr "${chr}" pfile_in
  pfile_name=chr${chr}_for_regenie_${stratum}

  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting regenie single variant conditioning task, version ${version}" > ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] phenotype_tsv: ${phenotype_tsv}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] covariate_tsv: ${covariate_tsv}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] exclcovars: ${exclcovars}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] conditional: ${conditional}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] type: ${type}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] invnorm: ${invnorm}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] catvars: ${catvars}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] chr: ${chr}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] vars_to_run: ${vars_to_run}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] condition_vars: ${condition_vars}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] stratum: ${stratum}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] pheno_to_run: ${pheno_to_run}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] check_collinearity: ${check_coll}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] cond_engine: ${cond_engine}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] run_prefix: ${run_prefix}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] step1_run_prefix: ${step1_run_prefix}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] run_suffix: ${run_suffix}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] step1_directory: ${step1_directory}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] step1_input_template: ${step1_input_template}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] results_directory: ${results_directory}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] pfile_template: ${pfile_template}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] stratify_covar: ${stratify_covar}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] stratify_value: ${stratify_value}" >> ${log_file}

  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Setting up environment" >> ${log_file}
  setup
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Checking if run has already completed" >> ${log_file}
  check_if_run
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Downloading input files" >> ${log_file}
  download
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Setting up conditional files" >> ${log_file}
  conditional_lists
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running regenie association testing" >> ${log_file}
  run_regenie
  if (( check_coll )); then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Testing for collinearity between test and cond vars" >> ${log_file}
    check_collinearity
  fi
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Uploading results" >> ${log_file}
  gcloud storage cp ${outputdir}/* ${bucket}/${results_directory}/${outputdir}/
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed regenie single variant conditioning task" >> ${log_file}
  gcloud storage cp ${log_file} ${bucket}/${results_directory}/${outputdir}/
}

setup() {
  # make output directory
  outputdir=$(wgs_output_base "${type}" "${invnorm}" "${run_prefix}")
  plinkpath=${outputdir}"/meta_analysis/${pheno_to_run}_singlevar_conditioning_${run_suffix}/plink_files"
  outputdir=${outputdir}"/meta_analysis/${pheno_to_run}_singlevar_conditioning_${run_suffix}/singles_cond_step${step}"
  mkdir -p ${outputdir}
}

check_if_run() {
  # we always upload log last, so can check for it to determine if run completed
  # if no log, then check for any output files and delete them
  if gcloud storage objects describe "${bucket}/${results_directory}/${outputdir}/${log_file}" >/dev/null 2>&1; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Log file ${log_file} already exists, assuming complete and exiting"
    exit 0
  else
    gcloud storage rm "${bucket}/${results_directory}/${outputdir}/${job_prefix}*" 2>/dev/null || true
    if [ $step -eq 0 ]; then
      gcloud storage rm "${bucket}/${results_directory}/${plinkpath}/${pfile_name}*" 2>/dev/null || true
    fi
  fi
}

download() {
  # download input files, unzip if required

  # download conditional list if required
  if [[ $conditional != "" ]]; then
    gcloud storage cp "$conditional" conditional_list
  else
    touch conditional_list
  fi

  gcloud storage cp "$vars_to_run" vars_to_run_unfiltered
  awk -v c=${chr} '$1==c {print $2}' vars_to_run_unfiltered > vars_to_run

  decompress_to "$covariate_tsv" covariate_tsv
  decompress_to "$phenotype_tsv" phenotype_tsv

  awk -F':' -v c=${chr} '$1==c {print $0}' conditional_list > conditional_list_chr
  cat additional_vars_conditioning >> conditional_list_chr

  gcloud storage cp \
    "${bucket}/${step1_directory}/${step1_run_prefix}_${type}_${invnorm}_Step1_${stratum}_*" \
    .
  fetch_step1_sample_file_gcs "${step1_input_directory}" "${step1_input_name}"
  
  # move regenie step 1 list file (source uses the step1 prefix, target uses this run's prefix)
  mv ${step1_run_prefix}_${type}_${invnorm}_Step1_${stratum}_pred.list ${run_prefix}_Step1_pred.list

  # process stratum IDs
  extract_stratum_ids "${step1_input_name}" > stratum_ids

  # subset covariate and phenotype files by stratum
  subset_to_ids stratum_ids covariate_tsv
  subset_to_ids stratum_ids phenotype_tsv
  # then by covariate value, if this is a covariate-stratified run
  subset_to_covar_value "${stratify_covar}" "${stratify_value}" covariate_tsv phenotype_tsv

  cat conditional_list_chr vars_to_run > inclusion_list_final
  if [ ${step} -eq 0 ]; then
    # download WGS data
    gcloud storage cp "${pfile_in_directory}/${pfile_in_name}.p*" .
    # subset to covariate and phenotype files
    build_subset_ids covariate_tsv phenotype_tsv > subset_ids
    set +e
    plink2 --pfile ${pfile_in_name} --mac 11 --keep subset_ids --extract inclusion_list_final --make-pgen --out ${pfile_name}
    rc=$?
    set -e
    if [[ $rc -ne 0 ]]; then
      if grep -q "Error: No variants remaining after main filters." ${pfile_name}.log; then
        echo NO_VARIANTS_REMAINED_AFTER_MAC_FILTER > ${job_prefix}_no_variants_remaining
        gcloud storage cp ${job_prefix}_no_variants_remaining ${bucket}/${results_directory}/${outputdir}/
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] PLINK reported no variants to condition above mac threshold in this stratum, exiting" >> ${log_file}
        gcloud storage cp ${log_file} ${bucket}/${results_directory}/${outputdir}/
        exit 0
      else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Unhandled error in plink2 command, exiting" >> ${log_file}
        exit $rc
      fi
    fi
    gcloud storage cp ${pfile_name}.* ${bucket}/${results_directory}/${plinkpath}/
    # add _subset_step suffix to pfile
    for extension in pgen psam pvar; do
      mv ${pfile_name}.${extension} ${pfile_name}_subset_step.${extension}
    done
  else
    # download pfile and subset to step, also download .ma for manc-cojo if present
    gcloud storage cp "${bucket}/${results_directory}/${plinkpath}/${pfile_name}.*" .
    if [ "${cond_engine}" == "mancojo" ] && [ ! -s ${pfile_name}.ma ]; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] No ${pfile_name}.ma in ${plinkpath}; re-run step 0 with cond_engine=mancojo" >> ${log_file}
      cat ${log_file} >&2
      exit 1
    fi
    awk 'FNR==NR{arr[$1];next} ($3 in arr){print $3}' inclusion_list_final ${pfile_name}.pvar > inclusion_list_subset
    n_inclusion_subset=$(wc -l < inclusion_list_subset)
    if (( n_inclusion_subset == 0 )); then
      echo NO_VARIANTS_REMAINED_AFTER_MAC_FILTER > ${job_prefix}_no_variants_remaining
      gcloud storage cp ${job_prefix}_no_variants_remaining ${bucket}/${results_directory}/${outputdir}/
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] PLINK reported no variants to condition above mac threshold in this stratum, exiting" >> ${log_file}
      gcloud storage cp ${log_file} ${bucket}/${results_directory}/${outputdir}/
      exit 0
    fi
    plink2 --pfile ${pfile_name} --extract inclusion_list_final --make-pgen --out ${pfile_name}_subset_step
  fi
}

conditional_lists() {
  # check how many variants to condition on
  awk 'FNR==NR{arr[$1];next} ($3 in arr){print $3}' conditional_list_chr ${pfile_name}_subset_step.pvar > conditional_list_local
  n_conditional_snps=$(wc -l < conditional_list_local)
  # Past step 0 we must be conditioning on at least one variant
  # but they may not have carriers in every stratum.
  # If no carriers, output unconditioned summary stats from regenie
  if (( step > 0 && n_conditional_snps == 0 )); then
    n_selected=$(wc -l < conditional_list_chr)
    echo "NO_CONDITIONING_VARIANTS_PRESENT_IN_STRATUM ${step} ${n_selected}" > ${job_prefix}_no_conditioning_variants
    gcloud storage cp ${job_prefix}_no_conditioning_variants ${bucket}/${results_directory}/${outputdir}/
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] None of the ${step} selected hit(s) are present in this stratum; running unconditioned" >> ${log_file}
  fi

  # check if any snps to test
  awk 'FNR==NR{arr[$1];next} ($3 in arr){print $3}' vars_to_run ${pfile_name}_subset_step.pvar > vars_to_run_in_pvar
  if (( n_conditional_snps == 0 )); then
   cp vars_to_run_in_pvar vars_to_run_in_pvar_excl_conditional
  else
    awk 'FNR==NR{arr[$1];next} !($1 in arr){print $1}' conditional_list_local vars_to_run_in_pvar > vars_to_run_in_pvar_excl_conditional
  fi
  n_snps_to_run=$(wc -l < vars_to_run_in_pvar_excl_conditional)

  # exit if no variants to test
  if (( n_snps_to_run == 0 )); then
    echo "NO_VARIANTS_REMAINED_AFTER_REMOVING_CONDITIONAL_VARIANTS" > ${job_prefix}_no_variants_remaining
    gcloud storage cp ${job_prefix}_no_variants_remaining ${bucket}/${results_directory}/${outputdir}/
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] No variants remained after removing conditional variants, exiting" >> ${log_file}
    gcloud storage cp ${log_file} ${bucket}/${results_directory}/${outputdir}/
    exit 0
  fi

  # check if there are too many conditional variants (for the cohort size)
  var_issue=0
  n_samples=$(wc -l < <(grep -vE '^#' ${pfile_name}_subset_step.psam))
  if (( n_conditional_snps >= 20000 )); then
    echo "TOO_MANY_VARIANTS_TO_CONDITION_ON ${n_conditional_snps} ${n_samples}" > ${job_prefix}_too_many_variants
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Too many variants to condition on (${n_conditional_snps}), exiting" >> ${log_file}
    var_issue=1
  elif (( n_conditional_snps > n_samples / 2 )); then
    echo "TOO_MANY_VARIANTS_TO_CONDITION_ON_FOR_SAMPLE_SIZE ${n_conditional_snps} ${n_samples}" > ${job_prefix}_too_many_variants
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Number of variants to condition on (${n_conditional_snps}) is over half of sample size (${n_samples}), exiting" >> ${log_file}
    var_issue=1
  elif (( n_conditional_snps > n_samples / 10 )); then
    echo "NUMBER_OF_CONDITIONAL_VARIANTS_OVER_ONE_TENTH_OF_SAMPLE_SIZE ${n_conditional_snps} ${n_samples}" > ${job_prefix}_many_variants
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Warning: number of variants to condition on (${n_conditional_snps}) is over one tenth of sample size (${n_samples}), results may be unreliable" >> ${log_file}
    gcloud storage cp ${job_prefix}_many_variants ${bucket}/${results_directory}/${outputdir}/
  fi
  if (( var_issue == 1 )); then
    gcloud storage cp ${job_prefix}_too_many_variants ${bucket}/${results_directory}/${outputdir}/
    gcloud storage cp ${log_file} ${bucket}/${results_directory}/${outputdir}/
    exit 0
  fi
}

run_regenie() {
  # From step 1 on, manc-cojo replaces regenie fit if chosen. Step 0 always
  # runs regenie: its results are the unconditioned baseline the .ma is built from.
  if [ "${cond_engine}" == "mancojo" ] && [ ${step} -gt 0 ]; then
    run_mancojo
    return
  fi

  # run regenie
  bash run_regenie_step2_single_variants.sh ${run_prefix} ${catvars} ${outputdir} ${chr} ${invnorm} ${n_conditional_snps} ${job_prefix} ${pfile_name}_subset_step None ${type} ${exclcovars} ${pheno_to_run} ${threads}

  # build .ma to upload alongside pgen for manc-cojo
  if [ "${cond_engine}" == "mancojo" ]; then
    build_step0_ma
  fi

  # copy conditional list for upload
  cp conditional_list_local ${outputdir}/${job_prefix}.conditional_list
}

# Build the .ma every later step conditions from, and save it next to the pgen.
#
# COJO projects *marginal* statistics through an LD matrix, so the .ma must come
# from an unconditional fit. The step-0 regenie run above is only unconditional
# when no external conditional list was supplied, otherwise it will have run
# with --condition-list and its output is already adjusted (and omits the
# variants that were conditioned on). In that case a second, unconditional pass is
# run purely to build the .ma.
build_step0_ma() {
  local ma_dir=${outputdir} ma_prefix=${job_prefix}
  if (( n_conditional_snps > 0 )); then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] External conditional list in play at step 0; running an extra unconditional pass for the manc-cojo .ma" >> ${log_file}
    ma_dir=mancojo_work
    ma_prefix=${job_prefix}_uncond
    rm -rf ${ma_dir}
    mkdir -p ${ma_dir}
    bash run_regenie_step2_single_variants.sh ${run_prefix} ${catvars} ${ma_dir} ${chr} ${invnorm} 0 \
      ${ma_prefix} ${pfile_name}_subset_step None ${type} ${exclcovars} ${pheno_to_run} ${threads}
  fi

  shopt -s nullglob
  local marginal_files=( ${ma_dir}/${ma_prefix}*.regenie.gz )
  shopt -u nullglob
  if (( ${#marginal_files[@]} == 0 )); then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] No regenie output to build the manc-cojo .ma from, exiting" >> ${log_file}
    cat ${log_file} >&2
    exit 1
  fi
  local n_ma
  n_ma=$(build_mancojo_ma ${pfile_name}.ma "${marginal_files[@]}")
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Built manc-cojo marginal input with ${n_ma} variants" >> ${log_file}
  gcloud storage cp ${pfile_name}.ma ${bucket}/${results_directory}/${plinkpath}/
}

run_mancojo() {
  # manc-cojo ignores conditioning IDs absent from the .ma, errors on an empty
  # list, and writes nothing when none match. With nothing to condition on, emit
  # the unconditional results instead, as the regenie engine does.
  local n_cond_in_ma=0
  if (( n_conditional_snps > 0 )); then
    n_cond_in_ma=$(awk 'NR==FNR{c[$1];next} FNR>1 && ($1 in c)' conditional_list_local ${pfile_name}.ma | wc -l)
  fi
  if (( n_cond_in_ma < n_conditional_snps )); then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $(( n_conditional_snps - n_cond_in_ma )) of ${n_conditional_snps} conditioning variants have no marginal row in ${pfile_name}.ma and will not be conditioned on" >> ${log_file}
  fi
  if (( n_cond_in_ma == 0 )); then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] No conditioning variants usable by manc-cojo; writing unconditional results" >> ${log_file}
    write_uncond_cma ${pfile_name}.ma ${pfile_name}_subset_step.pvar ${chr} ${outputdir}/${job_prefix}.cma.cojo
  else
    bash run_mancojo_conditional.sh \
      ${pfile_name}_subset_step \
      ${pfile_name}.ma \
      conditional_list_local \
      ${chr} \
      ${outputdir}/${job_prefix} \
      ${threads} \
      ${cond_rsq_threshold} \
      ${cond_mancojo_effect_mode}
  fi
  # manc-cojo writes nothing when every candidate is conditioned away or screened
  # out as collinear, so report that as no variants remaining
  if [[ ! -s ${outputdir}/${job_prefix}.cma.cojo ]]; then
    echo "NO_MANCOJO_CANDIDATES_REMAINING" > ${outputdir}/${job_prefix}_no_variants_remaining
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] manc-cojo returned no candidates, reporting no variants remaining" >> ${log_file}
  fi

  # copy conditional list for upload
  cp conditional_list_local ${outputdir}/${job_prefix}.conditional_list
}

check_collinearity() {
  # can only run if at least one cond var present in cohort (already exited if no test var)
  if (( n_conditional_snps == 0 )); then
    echo "NO_CONDITIONAL_VARIANTS" > ${outputdir}/${job_prefix}_no_collinearity
  else
    fedcat collinearity \
      --pfile ${pfile_name}_subset_step \
      --selected conditional_list_local \
      --candidate vars_to_run_in_pvar_excl_conditional \
      --out ${outputdir}/${job_prefix}_variant_collinearity \
      --threads ${threads}
  fi
}


main
