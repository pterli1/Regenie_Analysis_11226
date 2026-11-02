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
chrom_start=${9}
chrom_end=${10}
stratum=${11}
phenos_to_run=${12}
bucket=${13}
run_prefix=${14}
step1_directory=${15}
step1_input_template=${16}
results_directory=${17}
pfile_template=${18}
annotation_manifest=${19}
coding_only=${20}
aaf_bins=${21}
# Prefixes (comma-separated) under which step1 predictions were written; defaults to
# run_prefix so older callers that don't pass it keep working.
step1_run_prefix=${22:-$run_prefix}
# Covariate stratification: covariate column and the value to keep. Empty (the default
# for callers that don't pass them) means run on everyone.
stratify_covar=${23:-}
stratify_value=${24:-}
# Tag naming this phenotype group's outputs; empty (the default) for an ungrouped run.
pheno_group=${25:-}

source shared_functions.sh
source annotation_helpers.sh

handle_pipefails() {
  # ignore exit code 141 from more complex command pipes
  # - use with: cmd1 | cmd2 | cmd3 || handle_pipefails "${PIPESTATUS[@]}"
  # need this where head closes stdin so cat and sort etc. return 141 SIGPIPE
  for x in "$@"; do
    (( $x == 141 )) || { (( $x > 0 )) && return $x; }
  done
  return 0
}

main() {
  version="1.7.1"
  # resolve templates into <name>_path/_directory/_name
  split_template "${pfile_template}" chr "${chr}" pfile
  split_template "${step1_input_template}" stratum "${stratum}" step1_input
  # parse annotation_manifest into manifest_resolved TSV (resolving {chr})
  annots_running=$(parse_annotation_manifest "${annotation_manifest}" "${chr}")
  chunk_prefix=${run_prefix}_${stratum}_chr${chr}_${chrom_start}_${chrom_end}${pheno_group:+_${pheno_group}}_${annots_running}_step2_aggregates
  log_file=${chunk_prefix}.workflow.log
  threads=$(nproc)

  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting regenie step2 aggregates task, version ${version}" > ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] phenotype_tsv: ${phenotype_tsv}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] covariate_tsv: ${covariate_tsv}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] exclcovars: ${exclcovars}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] conditional: ${conditional}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] type: ${type}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] invnorm: ${invnorm}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] catvars: ${catvars}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] chr: ${chr}" >> ${log_file}
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
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] annotation_manifest: ${annotation_manifest}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] coding_only: ${coding_only}" >> ${log_file}
  
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Setting up environment" >> ${log_file}
  setup
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Checking if run has already completed" >> ${log_file}
  check_if_run
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Downloading input files" >> ${log_file}
  download
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Setting up masks" >> ${log_file}
  mask_setup
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Subsetting input files" >> ${log_file}
  subset_input_files
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running regenie association testing" >> ${log_file}
  run_regenie
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Uploading results" >> ${log_file}
  find ${outputdir} -maxdepth 1 -type f -printf "%f\n" > ${outputdir}/${chunk_prefix}.output_filelist
  gcloud storage cp ${outputdir}/* ${bucket}/${results_directory}/${outputdir}/
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed regenie step2 aggregates task" >> ${log_file}
  gcloud storage cp ${log_file} ${bucket}/${results_directory}/${outputdir}/
}

setup() {
  # make output directory
  outputdir=$(wgs_output_base "${type}" "${invnorm}" "${run_prefix}")"/meta_analysis/aggregates"
  mkdir -p $outputdir
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

  # build annotation download list from manifest (paths are full paths)
  while IFS=$'\t' read -r name mask_path annot_path setlist_path; do
    mkdir -p ./input_files/${name}
    gcloud storage cp \
      "${mask_path}" \
      "${annot_path}" \
      "${setlist_path}" \
      "./input_files/${name}/"
  done < manifest_resolved

  # build local-name manifest (basenames without .zst) for later use
  build_manifest_local

  mapfile -t step1_files < <(step1_sources "${bucket}/${step1_directory}" "${step1_run_prefix}" "${type}" "${invnorm}" "${stratum}")
  gcloud storage cp \
    "${pfile_directory}/${pfile_name}.p*" \
    "${step1_files[@]}" \
    .
  fetch_step1_sample_file_gcs "${step1_input_directory}" "${step1_input_name}"

  # decompress annotation/pfile .zst files in parallel
  shopt -s nullglob
  zst_files=( *.zst ./input_files/*/*.zst )
  shopt -u nullglob
  if ! decompress_zst_parallel "${zst_files[@]}"; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: One or more decompression tasks failed" >> ${log_file}
    exit 1
  fi

  # merge the pred lists of every step1 set holding this group's phenotypes
  merge_step1_pred_lists "${step1_run_prefix}" "${type}" "${invnorm}" "${stratum}" "${phenos_to_run}" ${run_prefix}_Step1_pred.list

  # process stratum IDs
  extract_stratum_ids "${step1_input_name}" > stratum_ids
}

mask_setup() {
  # for each annotation triplet in the manifest, filter its set_list by
  # the requested genomic range and collect annotation files
  n_sets=()
  > variants_to_include_masks
  while IFS=$'\t' read -r name mask annot setlist; do
    # check no dots in any mask names
    if cut -f1 "$mask" | grep -q '\.'; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: one or more invalid masks with '.' in name in masks file ${mask}" >> ${log_file}
      exit 1
    fi
    # filter set_list by range
    # optionally subset to only coding sets (bodge)
    if [[ "${coding_only}" == "true" ]]; then
      awk -v chr=${chr} -v start=${chrom_start} -v end=${chrom_end} \
        '$2==chr && $3>=start && $3<=end && $1 ~ /^coding_/' "$setlist" > "${setlist}_chunk"      
    else
      awk -v chr=${chr} -v start=${chrom_start} -v end=${chrom_end} \
        '$2==chr && $3>=start && $3<=end' "$setlist" > "${setlist}_chunk"
    fi
    # Subset annotations file based on sets, count number of sets in chunk
    awk 'FNR==NR {arr[$1]; next} $2 in arr' "${setlist}_chunk" "$annot" > "${annot}_chunk"
    count=$(cut -f2 "${annot}_chunk" | sort -u | wc -l)
    n_sets+=("$count")
    # extract variants to include in block
    cut -f1 "${annot}_chunk" >> variants_to_include_masks
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${count} sets for triplet ${name}: ${mask} / ${annot} / ${setlist}" >> ${log_file}
  done < manifest_local

  # exit with rc 0 if set count is 0 for all annot triplets
  for n in "${n_sets[@]}"; do (( n != 0 )) && break; done
  if (( n == 0 )) ; then
    touch "${chunk_prefix}.output_filelist"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] No sets in region, exiting (success)" >> ${log_file}
    gcloud storage cp "${chunk_prefix}.output_filelist" "${log_file}" "${bucket}/${results_directory}/${outputdir}/"
    exit 0
  fi
}

subset_input_files() {
  # subset conditional list to variants present on chr
  awk -F':' -v c=$chr '$1==c' conditional_list > conditional_list_chr
  if [ $(wc -l < conditional_list_chr) -eq 0 ]; then
    touch conditional_list_local
  else
    awk 'FNR==NR{arr[$1];next} ($3 in arr){print $3}' conditional_list_chr ${pfile_name}.pvar > conditional_list_local
  fi
  n_conditional_snps=$(wc -l < conditional_list_local)

  # subset covariate and phenotype files by stratum
  subset_to_ids stratum_ids covariate_tsv
  subset_to_ids stratum_ids phenotype_tsv
  # then by covariate value, if this is a covariate-stratified run
  subset_to_covar_value "${stratify_covar}" "${stratify_value}" covariate_tsv phenotype_tsv
  build_subset_ids covariate_tsv phenotype_tsv > subset_ids

  # make subset plink files
  if [[ ${n_conditional_snps} -gt 0 ]]; then
  	plink2 --pfile ${pfile_name} --extract conditional_list_local --keep subset_ids --make-pgen --double-id --out chr${chr}_for_conditioning --memory 25000
  fi

  # The block's annotated variants can all be absent from this stratum's pfile, or none
  # rare enough for the mac10 build; plink2 then stops with no variants left
  if ! plink2_allow_no_variants chr${chr}_for_regenie --pfile ${pfile_name} --extract variants_to_include_masks \
    --keep subset_ids --make-pgen --double-id --out chr${chr}_for_regenie --memory 25000; then
    touch "${chunk_prefix}.output_filelist"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] No annotated variants of this block in the stratum, exiting (success)" >> ${log_file}
    gcloud storage cp "${chunk_prefix}.output_filelist" "${log_file}" "${bucket}/${results_directory}/${outputdir}/"
    exit 0
  fi
  have_mac10=1
  plink2_allow_no_variants chr${chr}_for_regenie_mac10 --pfile ${pfile_name} --extract variants_to_include_masks \
    --keep subset_ids --max-mac 10 --make-pgen --double-id --out chr${chr}_for_regenie_mac10 --memory 25000 \
    || have_mac10=0
}

run_regenie() {
  # run regenie using xargs to run 2 jobs at once
  com=""
  sub_threads=$(( threads / 2 ))
  triplet_idx=0
  while IFS=$'\t' read -r name mask annot setlist; do
    if [ ${n_sets[$triplet_idx]} -gt 0 ]; then
      if (( have_mac10 )); then
        com="${com} $run_prefix $catvars $outputdir ${chunk_prefix}_${name} ${chr} $n_conditional_snps ${mask} ${annot}_chunk ${setlist}_chunk ${aaf_bins} $invnorm mac10 $type ${exclcovars} ${phenos_to_run} ${sub_threads}"
      fi
      com="${com} $run_prefix $catvars $outputdir ${chunk_prefix}_${name} ${chr} $n_conditional_snps ${mask} ${annot}_chunk ${setlist}_chunk ${aaf_bins} $invnorm full $type ${exclcovars} ${phenos_to_run} ${sub_threads}"
    fi
    triplet_idx=$((triplet_idx + 1))
  done < manifest_local
  echo ${com}
  echo $com | xargs -n 16 -P 2 bash -c 'bash run_regenie_step2_aggregates.sh "$@"' _

  # copy conditional list for upload
  cp conditional_list_local ${outputdir}/${chunk_prefix}.conditional_list
}


main
