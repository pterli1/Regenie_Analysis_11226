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
mac=${9}
chrom_start=${10}
chrom_end=${11}
stratum=${12}
phenos_to_run=${13}
bucket=${14}
run_prefix=${15}
step1_directory=${16}
step1_input_template=${17}
results_directory=${18}
pfile_template=${19}
# Prefixes (comma-separated) under which step1 predictions were written; defaults to
# run_prefix so older callers that don't pass it keep working.
step1_run_prefix=${20:-$run_prefix}
# Covariate stratification: covariate column and the value to keep. Empty (the default
# for callers that don't pass them) means run on everyone.
stratify_covar=${21:-}
stratify_value=${22:-}
# GxE interaction variable; empty (the default) means no interaction test.
interaction=${23:-}
# Tag naming this phenotype group's outputs; empty (the default) for an ungrouped run.
pheno_group=${24:-}

source shared_functions.sh

main() {
  version="1.7.1"
  # resolve templates into <name>_path/_directory/_name
  split_template "${pfile_template}" chr "${chr}" pfile
  split_template "${step1_input_template}" stratum "${stratum}" step1_input
  chunk_prefix=${run_prefix}_${stratum}_chr${chr}_${chrom_start}_${chrom_end}${pheno_group:+_${pheno_group}}_step2_singlevar
  log_file=${chunk_prefix}.workflow.log
  threads=$(nproc)

  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting regenie step2 single variants task, version ${version}" > ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] phenotype_tsv: ${phenotype_tsv}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] covariate_tsv: ${covariate_tsv}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] exclcovars: ${exclcovars}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] conditional: ${conditional}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] type: ${type}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] invnorm: ${invnorm}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] catvars: ${catvars}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] chr: ${chr}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] mac: ${mac}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] chrom_start: ${chrom_start}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] chrom_end: ${chrom_end}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] stratum: ${stratum}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] phenos_to_run: ${phenos_to_run}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] pheno_group: ${pheno_group:-none}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] run_prefix: ${run_prefix}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] step1_run_prefix: ${step1_run_prefix}" >> ${log_file}
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
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Subsetting input files" >> ${log_file}
  subset_input_files
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running regenie association testing" >> ${log_file}
  run_regenie
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Uploading results" >> ${log_file}
  find ${outputdir} -maxdepth 1 -type f -printf "%f\n" > ${outputdir}/${chunk_prefix}.output_filelist
  gcloud storage cp ${outputdir}/* ${bucket}/${results_directory}/${outputdir}/
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed regenie step2 single variants task" >> ${log_file}
  gcloud storage cp ${log_file} ${bucket}/${results_directory}/${outputdir}/
}

setup() {
  # make output directory
  outputdir=$(wgs_output_base "${type}" "${invnorm}" "${run_prefix}")"/meta_analysis/single_variants"
  mkdir -p ${outputdir}
}

check_if_run() {
  # we always upload log last, so can check for it to determine if run completed
  # if no log, then check for any output files and delete them
  if gcloud storage objects describe "${bucket}/${results_directory}/${outputdir}/${log_file}" >/dev/null 2>&1; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Log file ${log_file} already exists, assuming complete and exiting"
    exit 0
  else
    gcloud storage rm "${bucket}/${results_directory}/${outputdir}/${chunk_prefix}*" 2>/dev/null || true
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

  decompress_to "$covariate_tsv" covariate_tsv
  decompress_to "$phenotype_tsv" phenotype_tsv

  mapfile -t step1_files < <(step1_sources "${bucket}/${step1_directory}" "${step1_run_prefix}" "${type}" "${invnorm}" "${stratum}")
  gcloud storage cp \
    "${pfile_directory}/${pfile_name}.p*" \
    "${step1_files[@]}" \
    .
  fetch_step1_sample_file_gcs "${step1_input_directory}" "${step1_input_name}"

  # merge the pred lists of every step1 set holding this group's phenotypes
  merge_step1_pred_lists "${step1_run_prefix}" "${type}" "${invnorm}" "${stratum}" "${phenos_to_run}" ${run_prefix}_Step1_pred.list
  
  # process stratum IDs
  extract_stratum_ids "${step1_input_name}" > stratum_ids
}

subset_input_files() {
  # subset conditional list to variants present on chr
  awk -F':' -v c=$chr '$1==c' conditional_list > conditional_list_chr

  # generate list of variants to test based on position
  awk -v start_pos=${chrom_start} -v end_pos=${chrom_end} -v chr=${chr} '$2>=start_pos && $2<=end_pos {print $3}' ${pfile_name}.pvar > tmp_inclusion_list_temporary
  cat conditional_list_chr tmp_inclusion_list_temporary > inclusion_list_final

  # subset covariate and phenotype files by stratum
  subset_to_ids stratum_ids covariate_tsv
  subset_to_ids stratum_ids phenotype_tsv
  # then by covariate value, if this is a covariate-stratified run
  subset_to_covar_value "${stratify_covar}" "${stratify_value}" covariate_tsv phenotype_tsv
  build_subset_ids covariate_tsv phenotype_tsv > subset_ids

  # make a subset plink file
  plink2 --pfile ${pfile_name} --keep subset_ids --extract inclusion_list_final --mac ${mac} --make-pgen  --out chr${chr}_for_regenie --memory 25000
  
  # subset conditional list to variants present in pvar
  if [ $(wc -l < conditional_list_chr) -eq 0 ]; then
    touch conditional_list_local
  else
    awk 'FNR==NR{arr[$1];next} ($3 in arr){print $3}' conditional_list_chr chr${chr}_for_regenie.pvar > conditional_list_local
  fi
  n_conditional_snps=$(wc -l < conditional_list_local)
}

run_regenie() {
  # run regenie in series (conditional run takes longer so faster to run sequentially)
  bash run_regenie_step2_single_variants.sh ${run_prefix} ${catvars} ${outputdir} ${chr} ${invnorm} 0 ${chunk_prefix}_unconditional chr${chr}_for_regenie None ${type} ${exclcovars} ${phenos_to_run} ${threads} ${interaction}
  if [[ $conditional != "" ]]; then
    if [ ${n_conditional_snps} -ne 0 ]; then
      bash run_regenie_step2_single_variants.sh ${run_prefix} ${catvars} ${outputdir} ${chr} ${invnorm} ${n_conditional_snps} ${chunk_prefix} chr${chr}_for_regenie None ${type} ${exclcovars} ${phenos_to_run} ${threads} ${interaction}
    else
      # no variants to condition on for this chunk, so the conditional run would be
      # identical to the unconditional one — duplicate the unconditional output instead
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] No variants to condition on for chr${chr}:${chrom_start}-${chrom_end}, duplicating unconditional output for conditional results" >> ${log_file}
      uncond_prefix=${outputdir}/${chunk_prefix}_unconditional
      for f in "${uncond_prefix}"*; do
        suffix=${f#"${uncond_prefix}"}
        cp "$f" "${outputdir}/${chunk_prefix}${suffix}"
      done
    fi
  else
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] No conditional list provided for chr${chr}:${chrom_start}-${chrom_end}, skipping conditional run" >> ${log_file}
  fi

  # copy conditional list for upload
  cp conditional_list_local ${outputdir}/${chunk_prefix}.conditional_list
}

main
