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
aggregates=${9}
cond_aggregates=${10}
stratum=${11}
bucket=${12}
run_prefix=${13}
run_suffix=${14}
pheno_to_run=${15}
step1_directory=${16}
step1_input_template=${17}
results_directory=${18}
pfile_template=${19}
annotation_manifest=${20}
aaf_bins=${21}
# Prefix under which step1 predictions were written; defaults to run_prefix so
# older callers that don't pass it keep working.
step1_run_prefix=${22:-$run_prefix}
# Covariate stratification: covariate column and the value to keep. Empty (the default
# for callers that don't pass them) means run on everyone.
stratify_covar=${23:-}
stratify_value=${24:-}
# Binary traits only, and only for engines that fit a joint model: minimum number
# of minor allele copies required in each of cases and controls for a variant to
# be kept as a conditioning covariate. 1 drops only alleles confined to one arm,
# which separate the fit and are provably unidentifiable.
min_arm_ac=${25:-0}
# Engine used to derive conditional statistics: "regenie" re-fits an adjusted model
# per step; "mancojo" projects marginal statistics through an LD matrix.
cond_engine=${26:-regenie}
cond_rsq_threshold=${27:-0.9}
# Whether ACAT takes a stratum's single-variant result only where the variant is rare
# in it; off for callers that don't pass it, which keeps every result.
acat_singles_stratum_rare=${28:-false}
# manc-cojo --effect-mode (GCTA, imputeNA or removeNA)
cond_mancojo_effect_mode=${29:-imputeNA}
# Conditioning step number, from fedcat
cond_step=${30}
# Scratch for the unconditional regenie runs the manc-cojo engine needs; kept out
# of ${outputdir} so it is never uploaded or mistaken for this step's results.
mancojo_dir=mancojo_work

source shared_functions.sh
source annotation_helpers.sh

handle_pipefails() {
  # ignore exit code 141 from more complex command pipes
  # - use with: cmd1 | cmd2 | cmd3 || handle_pipefails "${PIPESTATUS[@]}"
  # need this where head closes stdin so cat and sort etc. return 141 SIGPIPE
  # Also allow exit code 1 for all the pesky grep's that match nothing!
  for x in "$@"; do
    (( x == 141 || x == 1 || x == 0)) || return $x;
  done
  return 0
}

main() {
  version="1.7.1"
  # resolve templates into <name>_path/_directory/_name
  split_template "${pfile_template}" chr "${chr}" pfile_in
  split_template "${step1_input_template}" stratum "${stratum}" step1_input
  pfile_name=chr${chr}_for_regenie_${stratum}
  # process aggregates to condition on
  if [ -z "$cond_aggregates" ]; then
    touch aggs_for_conditioning
  else
    echo "${cond_aggregates}" | sed 's/\,/\n/g' > aggs_for_conditioning
  fi
  # parse annotation_manifest into manifest_resolved TSV (resolving {chr})
  annots_running=$(parse_annotation_manifest "${annotation_manifest}" "${chr}")
  step=${cond_step}
  job_prefix="${pheno_to_run}_${stratum}_chr${chr}_${annots_running}_aggregate_conditioning_step${step}"
  log_file="${job_prefix}.workflow.log"
  threads=$(nproc)

  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting regenie aggregate conditioning task, version ${version}" > ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] phenotype_tsv: ${phenotype_tsv}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] covariate_tsv: ${covariate_tsv}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] exclcovars: ${exclcovars}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] conditional: ${conditional}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] type: ${type}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] invnorm: ${invnorm}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] catvars: ${catvars}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] chr: ${chr}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] aggregates: ${aggregates}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] cond_aggregates: ${cond_aggregates}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] stratum: ${stratum}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] pheno_to_run: ${pheno_to_run}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] run_prefix: ${run_prefix}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] step1_run_prefix: ${step1_run_prefix}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] run_suffix: ${run_suffix}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] step1_dir: ${step1_directory}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] step1_input_template: ${step1_input_template}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] results_dir: ${results_directory}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] pfile_template: ${pfile_template}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] stratify_covar: ${stratify_covar}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] stratify_value: ${stratify_value}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] annotation_manifest: ${annotation_manifest}" >> ${log_file}
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] cond_engine: ${cond_engine}" >> ${log_file}

  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Setting up environment" >> ${log_file}
  setup
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Checking if run has already completed" >> ${log_file}
  check_if_run
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Downloading input files" >> ${log_file}
  download
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Setting up masks" >> ${log_file}
  mask_setup
  expand_conditioning_aggregates
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Subsetting input files" >> ${log_file}
  subset_input_files
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running regenie association testing" >> ${log_file}
  aggregate_conditioning
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Uploading results" >> ${log_file}
  gcloud storage cp ${outputdir}/* ${bucket}/${results_directory}/${outputdir}/
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed regenie aggregate conditioning task" >> ${log_file}
  gcloud storage cp ${log_file} ${bucket}/${results_directory}/${outputdir}/
}

setup() {
  # make output directory
  outputdir=$(wgs_output_base "${type}" "${invnorm}" "${run_prefix}")
  plinkpath=${outputdir}"/meta_analysis/${pheno_to_run}_aggregate_conditioning_${run_suffix}/plink_files"
  outputdir=${outputdir}"/meta_analysis/${pheno_to_run}_aggregate_conditioning_${run_suffix}/agg_cond_step${step}"
  mkdir -p $outputdir
}

check_if_run() {
  # we always upload log last, so can check for it to determine if run completed
  # if no log, then check for any output files and delete them
  if gcloud storage objects describe "${bucket}/${results_directory}/${outputdir}/${log_file}" >/dev/null 2>&1; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Log file ${log_file} already exists, assuming complete and exiting"
    exit 0
  else
    gcloud storage rm "${bucket}/${results_directory}/${outputdir}/${job_prefix}*" 2>/dev/null || true
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

  gcloud storage cp "${aggregates}" aggs_to_test_unfiltered

  decompress_to "$covariate_tsv" covariate_tsv
  decompress_to "$phenotype_tsv" phenotype_tsv

  # download annotation triplets to subdirectories (paths are full paths)
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

  gcloud storage cp \
    "${bucket}/${step1_directory}/${step1_run_prefix}_${type}_${invnorm}_Step1_${stratum}_*" \
    .
  fetch_step1_sample_file_gcs "${step1_input_directory}" "${step1_input_name}"

  if [ ${step} -eq 0 ]; then
    gcloud storage cp "${pfile_in_directory}/${pfile_in_name}.p*" .
  else
    gcloud storage cp "${bucket}/${results_directory}/${plinkpath}/${pfile_name}.*" .
  fi

  # decompress annotation/pfile .zst files in parallel
  shopt -s nullglob
  zst_files=( *.zst ./input_files/*/*.zst )
  shopt -u nullglob
  if ! decompress_zst_parallel "${zst_files[@]}"; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: One or more decompression tasks failed" >> ${log_file}
    exit 1
  fi

  # process stratum IDs
  extract_stratum_ids "${step1_input_name}" > stratum_ids

  # move regenie step 1 list file (source uses the step1 prefix, target uses this run's prefix)
  mv ${step1_run_prefix}_${type}_${invnorm}_Step1_${stratum}_pred.list ${run_prefix}_Step1_pred.list

  # process aggs_to_test_unfiltered
  awk -v c=${chr} '$1==c {print $2}' aggs_to_test_unfiltered > assocs_for_current_round
  # extract SET names (everything before the last dot; masks cannot contain dots)
  sed 's/\.[^.]*$//' assocs_for_current_round > sets_for_current_round
  # also extract masks (everything after the last dot)
  sed 's/^.*\.//' assocs_for_current_round > masks_for_current_round
}

mask_setup() {
  # check, subset and extract variants from annot triplets
  > all_masks
  > variants_to_include_masks
  while IFS=$'\t' read -r name mask annot setlist; do
    # check no dots in any mask names
    if cut -f1 "$mask" | grep -q '\.'; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: one or more invalid masks with '.' in name in masks file ${mask}" >> ${log_file}
      exit 1
    fi
    # build all_masks file
    cat "$mask" >> all_masks
    # subset set list an annots to only those units in sets_for_current_round
    awk 'NR==FNR{arr[$1];next} $1 in arr' sets_for_current_round "$setlist" > "${setlist}_test_subset"
    awk 'NR==FNR{arr[$1];next} $2 in arr' sets_for_current_round "$annot" > "${annot}_test_subset"
    # extract variants in masks
    cut -f1 "${annot}_test_subset" >> variants_to_include_masks
  done < manifest_local

  # check that there are no duplicate mask names
  if [[ $(cut -f1 all_masks | sort | uniq -d) ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: duplicate mask names found across mask files" >> ${log_file}
    exit 1
  fi
}

subset_input_files() {
  # subset conditional list to variants present on chr
  awk -F':' -v c=$chr '$1==c' conditional_list > conditional_list_chr

  # subset covariate and phenotype files by stratum
  subset_to_ids stratum_ids covariate_tsv
  subset_to_ids stratum_ids phenotype_tsv
  # then by covariate value, if this is a covariate-stratified run
  subset_to_covar_value "${stratify_covar}" "${stratify_value}" covariate_tsv phenotype_tsv
  build_subset_ids covariate_tsv phenotype_tsv > subset_ids

  cat conditional_list_chr variants_to_include_masks condition_snps_temp > inclusion_list_final

  # create subset plink file for this step
  if [ ${step} -eq 0 ]; then
    # at step 0, subset from main input plink files
    set +e
    plink2 --pfile ${pfile_in_name} --mac 1 --keep subset_ids --extract inclusion_list_final --make-pgen --double-id --out ${pfile_name} --memory 25000
    rc=$?
    set -e
    if [[ $rc -ne 0 ]]; then
      if grep -q "Error: No variants remaining after main filters." ${pfile_name}.log; then
        echo "NO_POLYMORPHIC_VARIANTS_REMAINED_AFTER_FILTERING" > ${job_prefix}_no_aggregates_remaining
        gcloud storage cp ${job_prefix}_no_aggregates_remaining ${bucket}/${results_directory}/${outputdir}/
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] PLINK reported no polymorphic variants in this stratum, exiting" >> ${log_file}
        gcloud storage cp ${log_file} ${bucket}/${results_directory}/${outputdir}/
        exit 0
      else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Unhandled error in plink2 command, exiting" >> ${log_file}
        exit $rc
      fi
    fi
    gcloud storage cp ${pfile_name}.p* ${bucket}/${results_directory}/${plinkpath}/
    for extension in pgen psam pvar; do
      mv ${pfile_name}.${extension} chr${chr}_for_regenie.${extension}
    done
  else
    # cached pgen is already filtered to this stratum, filter varints to this step
    plink2 --pfile ${pfile_name} --extract inclusion_list_final --make-pgen --double-id --out chr${chr}_for_regenie --memory 25000
  fi

  # subset conditional list to variants present in stratum
  if [ $(wc -l < conditional_list_chr) -eq 0 ]; then
    touch conditional_list_local
  else
    awk 'FNR==NR{arr[$1];next} ($3 in arr){print $3}' conditional_list_chr chr${chr}_for_regenie.pvar > conditional_list_local
  fi
  n_conditional_snps=$(wc -l < conditional_list_local)

  # check whether there are mac10 burdens and/or single varints to run
  cat conditional_list_local condition_snps_temp > not_test_variants
  # Every remaining variant can be a conditioning variant (e.g. nested masks of a
  # selected aggregate), which would leave plink2 nothing to count
  n_test_variants=$(awk 'FILENAME == ARGV[1] {arr[$1]; next} !/^#/ && !($3 in arr)' not_test_variants chr${chr}_for_regenie.pvar | wc -l)
  if (( n_test_variants == 0 )); then
    echo "NO_AGGREGATES_REMAINED_AFTER_REMOVING_CONDITIONAL_VARIANTS" > ${job_prefix}_no_aggregates_remaining
    gcloud storage cp ${job_prefix}_no_aggregates_remaining ${bucket}/${results_directory}/${outputdir}/
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Every variant left in the stratum is a conditioning variant, exiting" >> ${log_file}
    gcloud storage cp ${log_file} ${bucket}/${results_directory}/${outputdir}/
    exit 0
  fi
  plink2 --pfile chr${chr}_for_regenie --exclude not_test_variants --freq counts cols=+reffreq,+alt1freq,-maybeprovref --nonfounders --out all_test_variants --memory 25000
  # col 5/6 REF_CT/ALT1_CT ref/alt counts
  n_mac_lteq_10=$(awk 'FNR > 1 && ($5 <= 10 || $6 <= 10)' all_test_variants.acount | wc -l)
  n_mac_gteq_11=$(awk 'FNR > 1 && ($5 >= 11 && $6 >= 11)' all_test_variants.acount | wc -l)

  # subset plink if required
  if (( n_mac_lteq_10 > 0 )); then
    plink2 --pfile chr${chr}_for_regenie --max-mac 10 --make-pgen --double-id --out chr${chr}_for_regenie_mac10 --memory 25000
  fi

  if (( n_mac_gteq_11 > 0 )); then
    plink2 --pfile chr${chr}_for_regenie --mac 11 --make-pgen --double-id --out chr${chr}_for_regenie_for_sv --memory 25000
  fi
}

# ── manc-cojo engine ─────────────────────────────────────────────────────────
# COJO derives conditional effects from marginal statistics, so the inputs it
# needs are the unconditional ones. Two things make aggregates harder than single
# variants here:
#
#  1. The units tested are burden masks, which are built inside regenie and are
#     not in any pgen. --write-mask materialises them so LD between a mask and the
#     variants being conditioned on can be computed.
#  2. The variants being conditioned on are an aggregate's rare constituents, and
#     regenie never reports them individually — they only ever appear inside a
#     mask. So a single-variant pass over every mask variant is needed purely to
#     give COJO their marginal effects.
#
# These are rebuilt at every step rather than cached. The mask definition itself
# changes as conditioning proceeds — ${annot}_test_subset_rare_nocond drops the
# variants being conditioned on, exactly as the regenie engine does — so a cached
# mask would be stale, and keeping the same definition is what makes the two
# engines answer the same question and be comparable.
build_mancojo_agg_inputs() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Building manc-cojo marginal inputs (unconditional)" >> ${log_file}

  # Scratch, NOT ${outputdir}: everything under ${outputdir} is uploaded to the
  # results directory and re-globbed by fedcat as this step's conditional
  # results. Unconditional regenie output matches those globs exactly, so
  # writing it there would let it be meta-analysed as if it were conditional.
  rm -rf ${mancojo_dir}
  mkdir -p ${mancojo_dir}/aggs ${mancojo_dir}/singles

  # 1. Unconditional aggregate runs, writing the mask genotypes out. ncon=0 makes
  #    these unconditional even when an external conditional list was supplied.
  local triplet_idx=0
  while IFS=$'\t' read -r name mask annot setlist; do
    if [ ${n_triplet_sets[$triplet_idx]:-0} -gt 0 ]; then
      if (( ${n_triplet_mac10_snps[$triplet_idx]:-0} > 0 )); then
        bash run_regenie_step2_aggregates.sh ${run_prefix} ${catvars} ${mancojo_dir}/aggs ${job_prefix}_uncond_${name} \
          ${chr} 0 ${mask} ${annot}_test_subset_rare_nocond ${setlist}_test_subset_aftercond ${aaf_bins} \
          ${invnorm} mac10 ${type} ${exclcovars} ${pheno_to_run} ${threads} 1
      fi
      bash run_regenie_step2_aggregates.sh ${run_prefix} ${catvars} ${mancojo_dir}/aggs ${job_prefix}_uncond_${name} \
        ${chr} 0 ${mask} ${annot}_test_subset_rare_nocond ${setlist}_test_subset_aftercond ${aaf_bins} \
        ${invnorm} full ${type} ${exclcovars} ${pheno_to_run} ${threads} 1
    fi
    triplet_idx=$((triplet_idx + 1))
  done < manifest_local

  # 2. Marginal statistics for every mask variant, including the ultra-rare ones
  #    that are only ever conditioned on. The trailing 1 is --minMAC: REGENIE
  #    defaults to 5 and the conditioning variants are overwhelmingly below that,
  #    and without a marginal row COJO cannot condition on them at all.
  bash run_regenie_step2_single_variants.sh ${run_prefix} ${catvars} ${mancojo_dir}/singles ${chr} ${invnorm} 0 \
    ${job_prefix}_uncond_allvars chr${chr}_for_regenie None ${type} ${exclcovars} ${pheno_to_run} ${threads} "" 1

  # 3. Assemble the .ma files: single-variant marginals plus each mask set.
  shopt -s nullglob
  local sv_files=( ${mancojo_dir}/singles/*.regenie.gz )
  local mac10_files=( ${mancojo_dir}/aggs/*_mac10_*.regenie.gz )
  local burden_files=( ${mancojo_dir}/aggs/*_burden_*.regenie.gz )
  # Mask filesets go via list files rather than arrays: an empty array expanded
  # under `set -u` is only safe from bash 4.4, and the task image is not pinned.
  # Empty filesets (no mask passed --minMAC) would fail plink2, so skip them.
  local bim
  : > mac10_mask_bims
  for bim in ${mancojo_dir}/aggs/*_mac10_masks.bim; do [[ -s "${bim}" ]] && echo "${bim}" >> mac10_mask_bims; done
  : > burden_mask_bims
  for bim in ${mancojo_dir}/aggs/*_burden_masks.bim; do [[ -s "${bim}" ]] && echo "${bim}" >> burden_mask_bims; done
  shopt -u nullglob

  if (( ${#sv_files[@]} == 0 )); then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] No unconditional single-variant output for the manc-cojo .ma, exiting" >> ${log_file}
    cat ${log_file} >&2
    exit 1
  fi

  build_mancojo_ma singles.ma "${sv_files[@]}" > /dev/null
  # A mask build can run and still test nothing: masks left below --minMAC once
  # the conditioning variants are removed give header-only output.
  local n_burden_rows=0 n_mac10_rows=0
  if (( ${#burden_files[@]} > 0 )); then
    n_burden_rows=$(build_mancojo_ma --allow-empty burden_masks.ma "${burden_files[@]}")
  fi
  if (( n_burden_rows > 0 )); then
    cat singles.ma > cojo_burden.ma
    tail -n +2 burden_masks.ma >> cojo_burden.ma
  fi
  if (( ${#mac10_files[@]} > 0 )); then
    n_mac10_rows=$(build_mancojo_ma --allow-empty mac10_masks.ma "${mac10_files[@]}")
  fi
  if (( n_mac10_rows > 0 )); then
    cat singles.ma > cojo_mac10.ma
    tail -n +2 mac10_masks.ma >> cojo_mac10.ma
  else
    # No mac10 masks to test. The mac10 build still runs, because it is the one
    # that carries the single-variant rows; the empty mask half is dropped below.
    cp singles.ma cojo_mac10.ma
  fi

  if (( n_burden_rows == 0 && n_mac10_rows == 0 )); then
    # Matches the regenie engine's header-only output: no aggregate is testable in
    # this stratum, but its single-variant results still go to the meta.
    echo "NO_MASKS_ABOVE_MINMAC_AFTER_CONDITIONING" > ${outputdir}/${job_prefix}_no_aggregates_remaining
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] No mask passed --minMAC after removing conditioning variants; writing single-variant results only" >> ${log_file}
  fi

  # manc-cojo's marginal inputs, kept for diagnostics; fedcat does not download them
  for build in mac10 burden; do
    if [[ -s cojo_${build}.ma ]]; then
      gzip -c cojo_${build}.ma > ${outputdir}/${job_prefix}_${build}_mancojo_input.ma.gz
    fi
  done

  # 4. Merge mask pseudo-genotypes into the variant pgen so manc-cojo can compute
  #    LD spanning both. Mask IDs collide between the mac10 and full builds (same
  #    SET.MASK.aafbin, different genotypes), so they are kept in separate filesets
  #    rather than merged into one.
  merge_mask_geno cojo_mac10_geno mac10_mask_bims
  if [[ -s cojo_burden.ma ]]; then
    merge_mask_geno cojo_burden_geno burden_mask_bims
  fi
}

# Merge the regenie mask filesets listed in a file into the mask-variant pgen.
#   merge_mask_geno <out prefix> <file listing mask .bim paths>
merge_mask_geno() {
  local out=$1 bim_list=$2
  if [ ! -s "${bim_list}" ]; then
    # no masks for this build; the variant pgen alone still supports conditioning
    plink2 --pfile chr${chr}_for_regenie --make-pgen --out ${out} --memory 25000
    return 0
  fi
  # plink2 --pmerge only handles filesets that follow one another in position
  # order, and masks interleave with the variants. Instead put every fileset on
  # the same samples in the same order, concatenate the .bed variant blocks
  # directly, then sort.
  sed 's/\.bim$/.fam/' "${bim_list}" | xargs awk '
    FNR == 1 { nf++ }
    { k = $1"\t"$2; c[k]++; if (nf == 1) o[++n] = k }
    END { for (i = 1; i <= n; i++) if (c[o[i]] == nf) print o[i] }' > ${out}_samples
  local inputs=( "--pfile chr${chr}_for_regenie" ) prefix i
  while read -r prefix; do inputs+=( "--bfile ${prefix%.bim}" ); done < "${bim_list}"
  : > ${out}_concat.bed
  : > ${out}_concat.bim
  for i in "${!inputs[@]}"; do
    # --indiv-sort needs pgen output, so reorder first and write .bed after
    plink2 ${inputs[$i]} --keep ${out}_samples --indiv-sort f ${out}_samples \
      --make-pgen --out ${out}_part${i}_tmp --memory 25000
    plink2 --pfile ${out}_part${i}_tmp --make-bed --out ${out}_part${i} --memory 25000
    if ! cmp -s <(cut -f1,2 ${out}_part0.fam) <(cut -f1,2 ${out}_part${i}.fam); then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] Sample order differs between ${out} mask filesets, exiting" >> ${log_file}
      cat ${log_file} >&2
      exit 1
    fi
    # keep the 3-byte .bed magic number from the first fileset only
    if (( i == 0 )); then cat ${out}_part${i}.bed; else tail -c +4 ${out}_part${i}.bed; fi >> ${out}_concat.bed
    cat ${out}_part${i}.bim >> ${out}_concat.bim
  done
  cp ${out}_part0.fam ${out}_concat.fam
  plink2 --bfile ${out}_concat --sort-vars --make-pgen --out ${out} --memory 25000
}

# One manc-cojo call per mask build, then split the result into the files fedcat
# already globs for: mask rows feed the burden/mac10 meta, variant rows the
# single-variant meta.
run_mancojo_aggregates() {
  build_mancojo_agg_inputs

  # The regenie engine only tests single variants with MAC >= 11, and only when
  # some exist. COJO returns an estimate for every non-conditioned row in the
  # .ma, so restrict it to the same set rather than emitting single-variant
  # results the regenie engine never would.
  : > sv_testable_ids
  if (( n_mac_gteq_11 > 0 )); then
    grep -vE '^#' chr${chr}_for_regenie_for_sv.pvar | cut -f3 | sort -u > sv_testable_ids
  fi

  local build bim n_cond_in_ma
  for build in mac10 burden; do
    [[ -s cojo_${build}.ma ]] || continue
    local cma=${outputdir}/${job_prefix}_${build}_mancojo.cma.cojo

    # manc-cojo silently ignores conditioning IDs absent from the .ma (regenie gave
    # no usable estimate), errors on an empty list, and writes nothing when none
    # match. With nothing to condition on, emit the unconditional results instead,
    # as the regenie engine does, so the stratum still reaches the meta.
    n_cond_in_ma=$(awk 'NR==FNR{c[$1];next} FNR>1 && ($1 in c)' conditional_list_local cojo_${build}.ma | wc -l)
    if (( n_cond_in_ma < n_conditional_snps )); then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] $(( n_conditional_snps - n_cond_in_ma )) of ${n_conditional_snps} conditioning variants have no marginal row in cojo_${build}.ma and will not be conditioned on" >> ${log_file}
    fi
    if (( n_cond_in_ma == 0 )); then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] No conditioning variants usable by manc-cojo for ${build}; writing unconditional results" >> ${log_file}
      write_uncond_cma cojo_${build}.ma cojo_${build}_geno.pvar ${chr} ${cma}
    else
      bash run_mancojo_conditional.sh \
        cojo_${build}_geno \
        cojo_${build}.ma \
        conditional_list_local \
        ${chr} \
        ${outputdir}/${job_prefix}_${build}_mancojo \
        ${threads} \
        ${cond_rsq_threshold} \
        ${cond_mancojo_effect_mode}
    fi
    [[ -s ${cma} ]] || continue

    # Mask rows are exactly the IDs in the mask bims; everything else is a single
    # variant. Splitting on the ID list avoids guessing from the ID format.
    : > mask_ids_${build}
    while IFS= read -r bim; do cut -f2 "${bim}"; done < ${build}_mask_bims >> mask_ids_${build}

    head -1 ${cma} > ${cma}.masks
    head -1 ${cma} > ${cma}.singles
    awk -v ids=mask_ids_${build} -v m=${cma}.masks -v v=${cma}.singles '
      BEGIN { while ((getline line < ids) > 0) is_mask[line] = 1 }
      FNR == 1 { next }
      { if ($2 in is_mask) print >> m; else print >> v }' ${cma}
    mv ${cma}.masks ${cma}
    # a header-only file is not "no result" to fedcat — it converts to an empty
    # frame that still enters the meta. Drop it so the build looks absent instead.
    if (( $(wc -l < ${cma}) <= 1 )); then
      rm -f ${cma}
    fi

    if [ "${build}" == "mac10" ] && (( n_mac_gteq_11 > 0 )); then
      # single-variant rows are identical across builds; keep the mac10 copy only
      head -1 ${cma}.singles > singles_filtered
      awk -v max="${singles_max_a1freq}" '
        NR == FNR { keep[$1]; next }
        FNR == 1 {
          for (i = 1; i <= NF; i++) if ($i == "freq") fc = i
          if (max != "" && !fc) { print "FATAL: no freq column in manc-cojo output" > "/dev/stderr"; exit 1 }
          next
        }
        ($2 in keep) && (max == "" || $fc + 0 <= max + 0)' \
        sv_testable_ids ${cma}.singles >> singles_filtered
      if (( $(wc -l < singles_filtered) > 1 )); then
        mv singles_filtered ${outputdir}/${job_prefix}_singlevars_mancojo.cma.cojo
      fi
    fi
    rm -f ${cma}.singles
  done

  # manc-cojo writes nothing when every candidate is conditioned away or screened
  # out as collinear; with no result and no marker fedcat would see a failed stratum
  if ! compgen -G "${outputdir}/${job_prefix}_*mancojo.cma.cojo" > /dev/null \
    && [[ ! -e ${outputdir}/${job_prefix}_no_aggregates_remaining ]]; then
    echo "NO_MANCOJO_CANDIDATES_REMAINING" > ${outputdir}/${job_prefix}_no_aggregates_remaining
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] manc-cojo returned no candidates, reporting no aggregates remaining" >> ${log_file}
  fi

  cp conditional_list_local ${outputdir}/${job_prefix}.conditional_list
}

# Keep annotation rows whose variant is present and rare in this stratum and is
# not being conditioned on.
filter_annots_rare_nocond() {
  local annots_in=$1 annots_out=$2 rc
  cut -f1 "${annots_in}" | sort -u > reformat_vars
  : > reformat_plink.snplist
  if [ -s reformat_vars ]; then
    set +e
    plink2 --pfile chr${chr}_for_regenie --extract reformat_vars --max-maf ${rare_maf} \
      --write-snplist --out reformat_plink --memory 25000
    rc=$?
    set -e
    # none present and rare in this stratum is an empty result, not a failure
    if (( rc != 0 )) && ! grep -q "No variants remaining" reformat_plink.log; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] Unhandled error in plink2 command, exiting" >> ${log_file}
      exit $rc
    fi
  fi
  awk 'NR==FNR{arr[$1];next} $1 in arr' reformat_plink.snplist "${annots_in}" > annotations_filtered_rare
  if [ ! -s all_condition_snps_unfiltered ]; then
    cp annotations_filtered_rare "${annots_out}"
  else
    awk 'NR==FNR{arr[$1];next} !($1 in arr)' all_condition_snps_unfiltered annotations_filtered_rare > "${annots_out}"
  fi
}

expand_conditioning_aggregates() {
  # Build list of variants to condition on from the aggregates being conditioned on.
  # Each entry in aggs_for_conditioning is SET.MASK — we split on the last dot.
  > condition_snps_temp
  while IFS=$'\n' read -r aggregate_temp; do
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Adding ${aggregate_temp} to adjustment list" >> ${log_file}
    top_mask="${aggregate_temp##*.}"
    top_set="${aggregate_temp%.*}"
    # masks are unique across triplets, so check each for match and only pull out vars from matching annots
    while IFS=$'\t' read -r name mask annot setlist; do
      awk -v mask="${top_mask}" '$1==mask {print $2}' "$mask" | tr ',' '\n' > mask_consq_temp
      n_mask_consq=$(wc -l < mask_consq_temp)
      if (( n_mask_consq > 0 )); then
        awk -v top_set="${top_set}" 'NR==FNR{arr[$1];next} $2==top_set && $3 in arr {print $1}' mask_consq_temp "$annot" | sort -u >> condition_snps_temp
        break
      fi
    done < manifest_local
  done < aggs_for_conditioning
}

aggregate_conditioning() {
  cat conditional_list_local condition_snps_temp > all_condition_snps_unfiltered

  # Mask variants are capped at 5x the largest AAF bin
  rare_maf=$(echo "${aaf_bins}" | tr ',' '\n' | awk 'NR == 1 || $1 > m { m = $1 } END { m *= 5; print (m > 0.5 ? 0.5 : m) }')
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Mask variant MAF cap: ${rare_maf}" >> ${log_file}

  # Single variants here only feed ACAT. When it takes a stratum's result only for
  # bins the variant is rare in, anything above the widest bin can never be used;
  # dropping it lets a stratum with nothing usable return no results (and stop
  # conditioning) instead of running on. Empty keeps every result.
  singles_max_a1freq=""
  if [[ "${acat_singles_stratum_rare}" == "true" ]]; then
    singles_max_a1freq=$(echo "${aaf_bins}" | tr ',' '\n' | awk 'NR == 1 || $1 > m { m = $1 } END { print m }')
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Single-variant results kept up to A1FREQ ${singles_max_a1freq}" >> ${log_file}
  fi

  # Check which aggs we are testing still have variants after removing conditional variants
  n_triplet_sets=()
  n_triplet_mac10_snps=()
  total_sets=0
  triplet_idx=0
  while IFS=$'\t' read -r name mask annot setlist; do
    # Only need to make extra checks when few aggs for the current triplet
    n_check=$(awk 'NR==FNR{arr[$1];next} $1 in arr' "$mask" masks_for_current_round | wc -l)
    if (( n_check == 0 )); then
      n_triplet_sets[$triplet_idx]=0
      n_triplet_mac10_snps[$triplet_idx]=0
    elif (( n_check > 30 )); then
      # if many aggs for the triplet, assume they will all have variants remaining and skip checks to save time
      n_triplet_sets[$triplet_idx]=$n_check
      n_triplet_mac10_snps[$triplet_idx]=999999
      total_sets=$((total_sets + n_triplet_sets[$triplet_idx]))
      filter_annots_rare_nocond "${annot}_test_subset" "${annot}_test_subset_rare_nocond"
      mv "${setlist}_test_subset" "${setlist}_test_subset_aftercond"
    else
      # For each SET.MASK in assocs_for_current_round, expand to SET<tab>CONSEQUENCE
      # pairs using the triplets mask definition file, so we can match against annotations
      # (which are VARIANT<tab>SET<tab>CONSEQUENCE).
      > assocs_for_current_triplet_with_consq
      while IFS=$'\n' read -r assoc; do
        set_temp="${assoc%.*}"
        mask_temp="${assoc##*.}"
        awk -v mask="${mask_temp}" '$1==mask {print $2}' "$mask" | tr ',' '\n' > consq_file
        while IFS=$'\n' read -r conseq; do
          printf '%s\t%s\n' "${set_temp}" "${conseq}" >> assocs_for_current_triplet_with_consq
        done < consq_file
      done < assocs_for_current_round

      # Subset annots to only rows matching a SET+CONSEQUENCE pair we care about.
      awk -F'\t' 'NR==FNR{pairs[$1 FS $2]; next} ($2 FS $3) in pairs' \
        assocs_for_current_triplet_with_consq "${annot}_test_subset" > annotations_filtered \
        || { code=$?; (( code == 1 )) || return $code; }

      filter_annots_rare_nocond annotations_filtered "${annot}_test_subset_rare_nocond"
      # work out which SETs still have variants remaining
      cut -f2 "${annot}_test_subset_rare_nocond" | sort -u > sets_with_variants_remaining
      # count how many of the assocs to test still have variants
      n_triplet_sets[$triplet_idx]=$(awk 'NR==FNR{arr[$1];next} ($1 in arr)' sets_with_variants_remaining assocs_for_current_triplet_with_consq | wc -l)
      total_sets=$((total_sets + n_triplet_sets[$triplet_idx]))

      # further subset the setlist for regenie
      awk 'NR==FNR{arr[$1];next} $1 in arr' sets_with_variants_remaining "${setlist}_test_subset" > "${setlist}_test_subset_aftercond"

      # check if any mac10 aggregates left as well
      if (( n_triplet_sets[$triplet_idx] > 0 )) && (( n_mac_lteq_10 > 0 )); then
        awk 'FNR==NR{arr[$1];next} ($3 in arr){print $3}' \
          "${annot}_test_subset_rare_nocond" \
          chr${chr}_for_regenie_mac10.pvar > triplet_mac10_vars_nocond
        n_triplet_mac10_snps[$triplet_idx]=$(wc -l < triplet_mac10_vars_nocond)
      else 
        n_triplet_mac10_snps[$triplet_idx]=0
      fi
    fi
    triplet_idx=$((triplet_idx + 1))
  done < manifest_local

  if (( total_sets == 0 )); then
    echo "NO_AGGREGATES_REMAINED_AFTER_REMOVING_CONDITIONAL_VARIANTS" > ${job_prefix}_no_aggregates_remaining
    gcloud storage cp ${job_prefix}_no_aggregates_remaining ${bucket}/${results_directory}/${outputdir}/
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] No aggregates remained after removing conditional variants, exiting" >> ${log_file}
    gcloud storage cp ${log_file} ${bucket}/${results_directory}/${outputdir}/
    exit 0
  fi

  # re-calulate number of conditional variants, including those from aggs
  if [ $(wc -l < all_condition_snps_unfiltered) -eq 0 ]; then
    n_conditional_snps=0
    cp conditional_list_local conditional_list_local_unfiltered
  else
    awk 'FNR==NR{arr[$1];next} ($3 in arr){print $3}' all_condition_snps_unfiltered chr${chr}_for_regenie.pvar > conditional_list_local
    # Any conditioning variant polymorphic in this stratum but missing from the step
    # pgen would be silently left out of the conditioning set.
    if (( step > 0 )); then
      awk 'FNR==NR{arr[$1];next} ($3 in arr){print $3}' all_condition_snps_unfiltered ${pfile_name}.pvar \
        | sort -u | comm -23 - <(sort -u conditional_list_local) > conditional_vars_missing
      if [ -s conditional_vars_missing ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Error: $(wc -l < conditional_vars_missing) conditioning variant(s) polymorphic in this stratum are missing from the step pgen, exiting" >> ${log_file}
        cat ${log_file} >&2
        exit 1
      fi
    fi
    # Pre-filter copy. Everything selected for conditioning must stay out of the
    # single-variant test set below, whether or not it survives as a covariate:
    # ACAT maps single-variant results back onto sets through the full step0
    # annotations, so a dropped variant that got tested would carry its marginal
    # p-value straight into the ACAT of the aggregate being conditioned on.
    cp conditional_list_local conditional_list_local_unfiltered
    # For binary traits, expanding an aggregate into its constituent variants puts
    # hundreds of ultra-rare covariates into a logistic model. An allele confined
    # to one arm separates the fit without contributing any case/control
    # information, so drop those as covariates. They stay excluded from testing
    # regardless: the single-variant set below and the annotation subsetting above
    # both key off the pre-filter list.
    #
    # Opt-in: the filter is off at the default of 0 and only runs once a positive
    # threshold is set. Only engines that fit a joint model need it at all -- COJO
    # projects marginal estimates through an LD matrix, so separation cannot arise
    # and the filter would only discard covariates it adjusts for perfectly well.
    if [ "${type}" == "Binary" ] && [ "${cond_engine}" != "mancojo" ] && (( min_arm_ac > 0 )); then
      filter_separating_conditional_vars \
        "chr${chr}_for_regenie" phenotype_tsv "${pheno_to_run}" "${min_arm_ac}" \
        conditional_list_local conditional_list_local_kept subset_ids
      n_dropped_separating=$(wc -l < conditional_list_local_kept.dropped)
      if (( n_dropped_separating > 0 )); then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Dropped ${n_dropped_separating} conditioning variants with fewer than ${min_arm_ac} minor allele copies in cases or in controls" >> ${log_file}
        cp conditional_list_local_kept.dropped ${outputdir}/${job_prefix}.conditional_list_dropped
      fi
      mv conditional_list_local_kept conditional_list_local
    elif [ "${type}" == "Binary" ] && [ "${cond_engine}" == "mancojo" ]; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] cond_engine=${cond_engine} fits no joint model; skipping the separating-covariate filter" >> ${log_file}
    elif [ "${type}" == "Binary" ]; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] min_arm_ac=0; separating-covariate filter disabled" >> ${log_file}
    fi
    n_conditional_snps=$(wc -l < conditional_list_local)
  fi
  
  # Past step 0 we must be conditioning on at least one aggregate
  # but they may not have carriers in every stratum.
  # If no carriers, output unconditioned summary stats from regenie
  if (( step > 0 && n_conditional_snps == 0 )); then
    n_before_filter=$(wc -l < conditional_list_local_unfiltered)
    if (( n_before_filter > 0 )); then
      # Present in the stratum but removed by min ac case/control filtering.
      echo "ALL_CONDITIONING_VARIANTS_DROPPED_BY_ARM_AC_FILTER ${step} ${n_before_filter}" > ${job_prefix}_no_conditioning_variants
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] All ${n_before_filter} conditioning variants from ${step} selected aggregate(s) were dropped by the min_arm_ac filter; running unconditioned" >> ${log_file}
    else
      n_expanded=$(wc -l < all_condition_snps_unfiltered)
      echo "NO_CONDITIONING_VARIANTS_PRESENT_IN_STRATUM ${step} ${n_expanded}" > ${job_prefix}_no_conditioning_variants
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] None of the ${n_expanded} variant(s) from ${step} selected aggregate(s) are polymorphic in this stratum; running unconditioned" >> ${log_file}
    fi
    gcloud storage cp ${job_prefix}_no_conditioning_variants ${bucket}/${results_directory}/${outputdir}/
  fi

  # Only need case count checks below for regenie, at least for testing
  if [ "${cond_engine}" == "mancojo" ] && [ ${step} -gt 0 ]; then
    run_mancojo_aggregates
    return
  fi

	#run sig hits conditional on top hit
	if [[ ${n_conditional_snps} -gt 0 ]]; then
    plink2 --pfile chr${chr}_for_regenie --extract conditional_list_local --make-pgen --double-id --out chr${chr}_for_conditioning --memory 25000
  fi

  # Check if there are too many conditional variants (for the cohort size)
  var_issue=0
  n_samples=$(wc -l < <(grep -vE '^#' chr${chr}_for_regenie.psam))
  # For a binary trait the limiting quantity is the number of cases, not the
  # number of samples: a conditioning set that looks small next to N can still be
  # large next to the case count, and it is the case count that drives separation.
  n_cases=0
  if [ "${type}" == "Binary" ]; then
    n_cases=$(comm -12 \
      <(build_case_ids phenotype_tsv "${pheno_to_run}" | cut -f2 | sort -u) \
      <(extract_stratum_ids "chr${chr}_for_regenie" | sort -u) | wc -l)
  fi
  if (( n_conditional_snps >= 20000 )); then
    echo "TOO_MANY_VARIANTS_TO_CONDITION_ON ${n_conditional_snps} ${n_samples}" > ${job_prefix}_too_many_variants
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Too many variants to condition on (${n_conditional_snps}), exiting" >> ${log_file}
    var_issue=1
  elif (( n_conditional_snps > n_samples / 2 )); then
    echo "TOO_MANY_VARIANTS_TO_CONDITION_ON_FOR_SAMPLE_SIZE ${n_conditional_snps} ${n_samples}" > ${job_prefix}_too_many_variants
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Number of variants to condition on (${n_conditional_snps}) is over half of sample size (${n_samples}), exiting" >> ${log_file}
    var_issue=1
  elif (( n_cases > 0 && n_conditional_snps > n_cases / 10 )); then
    echo "TOO_MANY_VARIANTS_TO_CONDITION_ON_FOR_CASE_COUNT ${n_conditional_snps} ${n_cases}" > ${job_prefix}_too_many_variants
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Number of variants to condition on (${n_conditional_snps}) is over one tenth of the case count (${n_cases}), exiting" >> ${log_file}
    var_issue=1
  elif (( n_cases > 0 && n_conditional_snps > n_cases / 20 )); then
    echo "NUMBER_OF_CONDITIONAL_VARIANTS_OVER_ONE_TWENTIETH_OF_CASE_COUNT ${n_conditional_snps} ${n_cases}" > ${job_prefix}_many_variants
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Warning: number of variants to condition on (${n_conditional_snps}) is over one twentieth of the case count (${n_cases}), results may be unreliable" >> ${log_file}
    gcloud storage cp ${job_prefix}_many_variants ${bucket}/${results_directory}/${outputdir}/
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

	# run regenie using xargs to run 2 jobs at once
  com=""
  sub_threads=$(( threads / 2 ))
  triplet_idx=0
  while IFS=$'\t' read -r name mask annot setlist; do
    if [ ${n_triplet_sets[$triplet_idx]} -gt 0 ]; then
      if [ ${n_triplet_mac10_snps[$triplet_idx]} -gt 0 ]; then
        com="${com} $run_prefix $catvars $outputdir ${job_prefix}_${name} ${chr} $n_conditional_snps ${mask} ${annot}_test_subset_rare_nocond ${setlist}_test_subset_aftercond ${aaf_bins} $invnorm mac10 $type ${exclcovars} ${pheno_to_run} ${sub_threads}"
      fi
      com="${com} $run_prefix $catvars $outputdir ${job_prefix}_${name} ${chr} $n_conditional_snps ${mask} ${annot}_test_subset_rare_nocond ${setlist}_test_subset_aftercond ${aaf_bins} $invnorm full $type ${exclcovars} ${pheno_to_run} ${sub_threads}"
    fi
    triplet_idx=$((triplet_idx + 1))
  done < manifest_local
  echo ${com}
  echo $com | xargs -n 16 -P 2 bash -c 'bash run_regenie_step2_aggregates.sh "$@"' _

  if (( n_mac_gteq_11 > 0 )); then
    # check if there are still single variants left to test after removing conditional variants
    if [ $(wc -l < conditional_list_local_unfiltered) -eq 0 ]; then
      grep -vE '^#' chr${chr}_for_regenie_for_sv.pvar | cut -f3 > single_vars_not_conditional
    else
      awk 'FNR==NR{arr[$1];next} !($3 in arr){print $3}' conditional_list_local_unfiltered <(grep -vE '^#' chr${chr}_for_regenie_for_sv.pvar) > single_vars_not_conditional
    fi
    # cond_pgen should be "None" no variants to condition on
    if (( n_conditional_snps > 0 )); then
      cond_pgen="chr${chr}_for_conditioning"
    else
      cond_pgen="None"
    fi
    n_single_vars=$(wc -l < single_vars_not_conditional)
    if (( n_single_vars > 0 )); then
      bash run_regenie_step2_single_variants.sh ${run_prefix} ${catvars} ${outputdir} ${chr} ${invnorm} ${n_conditional_snps} ${job_prefix}_singlevars chr${chr}_for_regenie_for_sv ${cond_pgen} ${type} ${exclcovars} ${pheno_to_run} ${threads}
      if [[ -n "${singles_max_a1freq}" ]]; then
        local sv_file n_kept
        for sv_file in ${outputdir}/${job_prefix}_singlevars_*.regenie.gz; do
          [[ -e "${sv_file}" ]] || continue
          n_kept=$(filter_regenie_max_a1freq "${sv_file}" "${singles_max_a1freq}")
          echo "[$(date '+%Y-%m-%d %H:%M:%S')] Kept ${n_kept} single-variant results with A1FREQ <= ${singles_max_a1freq}" >> ${log_file}
        done
      fi
    fi
  fi

  # copy conditional list for upload
  cp conditional_list_local ${outputdir}/${job_prefix}.conditional_list
}


main
