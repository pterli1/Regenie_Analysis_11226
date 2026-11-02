version 1.0

workflow meta_analysis_discovery_aou {

  input {
    Array[String] input_folders
    Array[String] input_biobanks
    Array[String] stratum_biobank_pairs
    Array[String] chromosomes
    Array[String] phenotypes
    String output_folder
    String run_prefix
    # Annotation source — provide exactly ONE of these two:
    #   schema_directory  → full GCS path to pre-built optimised parquet schema
    #   annotation_manifest → semicolon-separated AnnotationEntry triplets:
    #       name,mask_file,annot_file,setlist_file;name2,...
    #       Paths should be full GCS paths with {chr} placeholders.
    String schema_directory = ""
    String annotation_manifest = ""
    String dependencies_path = "gs://fc-secure-ab157c8b-aa91-4d13-b17f-8276c7b7fbcb/workflows/regenie_wgs_workflow/v1.7.1/dependencies"
    Boolean create_step0_subset = true
    Boolean run_uncond_singles = true
    Boolean run_cond_singles = true
    Boolean run_burden = true
    Boolean run_acat = true
    Boolean aggs_uncond = false
    Boolean acat_full_bounds = true
    String aaf_bins = "0.001"
    Boolean acat_singles_stratum_rare = true
    Float cond_p_singles_retest = 0.000001
    Float cond_p_agg_retest = 0.000001
    # Where each type of step-0 subset goes, with a completion marker per chromosome (default: output_folder)
    String step0_singles_folder = ""
    String step0_agg_folder = ""
    Int n_cpu = 16
    Int mem_gb = 128
    String disk = "200 SSD"
    Int max_retries = 0
    Int preemptible_tries = 0
}

  scatter (chr in chromosomes) {
    call meta_per_chromosome {
      input:
        chr = chr,
        input_folders = input_folders,
        input_biobanks = input_biobanks,
        stratum_biobank_pairs = stratum_biobank_pairs,
        phenotypes = phenotypes,
        output_folder = output_folder,
        run_prefix = run_prefix,
        schema_directory = schema_directory,
        annotation_manifest = annotation_manifest,
        dependencies_path = dependencies_path,
        create_step0_subset = create_step0_subset,
        run_uncond_singles = run_uncond_singles,
        run_cond_singles = run_cond_singles,
        run_burden = run_burden,
        run_acat = run_acat,
        aggs_uncond = aggs_uncond,
        acat_full_bounds = acat_full_bounds,
        aaf_bins = aaf_bins,
        acat_singles_stratum_rare = acat_singles_stratum_rare,
        cond_p_singles_retest = cond_p_singles_retest,
        cond_p_agg_retest = cond_p_agg_retest,
        step0_singles_folder = step0_singles_folder,
        step0_agg_folder = step0_agg_folder,
        n_cpu = n_cpu,
        mem_gb = mem_gb,
        disk = disk,
        max_retries = max_retries,
        preemptible_tries = preemptible_tries
    }
  }
}

task meta_per_chromosome {
  input {
    String chr
    Array[String] input_folders
    Array[String] input_biobanks
    Array[String] stratum_biobank_pairs
    String output_folder
    String run_prefix
    String schema_directory
    String annotation_manifest
    String dependencies_path
    Array[String] phenotypes
    Boolean create_step0_subset
    Boolean run_uncond_singles
    Boolean run_cond_singles
    Boolean run_burden
    Boolean run_acat
    Boolean aggs_uncond
    Boolean acat_full_bounds
    String aaf_bins
    Boolean acat_singles_stratum_rare
    Float cond_p_singles_retest
    Float cond_p_agg_retest
    String step0_singles_folder
    String step0_agg_folder
    Int n_cpu
    Int mem_gb
    String disk
    Int max_retries
    Int preemptible_tries
}

  command <<<
    set -euxo pipefail
    version="1.7.1"
    log_file=~{run_prefix}_meta_chr~{chr}.log
    gcloud config set storage/parallel_composite_upload_enabled False
    gcloud storage cp ~{dependencies_path}/* .
    chmod 777 *

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting meta-analysis task, version ${version}" > ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] chr: ~{chr}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] input_folders: ~{sep="," input_folders}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] input_biobanks: ~{sep="," input_biobanks}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] stratum_biobank_pairs: ~{sep="," stratum_biobank_pairs}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] output_folder: ~{output_folder}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] run_prefix: ~{run_prefix}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] schema_directory: ~{schema_directory}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] annotation_manifest: ~{annotation_manifest}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] phenotypes: ~{sep="," phenotypes}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] create_step0_subset: ~{create_step0_subset}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] cond_p_singles_retest: ~{cond_p_singles_retest}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] cond_p_agg_retest: ~{cond_p_agg_retest}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] step0_singles_folder: ~{step0_singles_folder}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] step0_agg_folder: ~{step0_agg_folder}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] n_cpu: ~{n_cpu}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] mem_gb: ~{mem_gb}" >> ${log_file}

    # Skip if this chromosome already completed. We upload the log last, so its
    # presence in the output folder means the run finished (independent of caching).
    if gcloud storage objects describe "~{output_folder}/${log_file}" >/dev/null 2>&1; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] Log file ${log_file} already exists, assuming complete and exiting"
      exit 0
    fi

    input_folders=(~{sep=' ' input_folders})
    input_biobanks=(~{sep=' ' input_biobanks})
    phenotypes=(~{sep=' ' phenotypes})

    do_uncond_singles="~{run_uncond_singles}"
    do_cond_singles="~{run_cond_singles}"
    do_burden="~{run_burden}"
    do_acat="~{run_acat}"
    aggs_uncond="~{aggs_uncond}"
    do_mac10=$do_acat

    # download discovery results for this chromosome.
    # fedcat reads .regenie.gz directly — no decompression or reformatting needed.
    for i in "${!input_biobanks[@]}"; do
      bb="${input_biobanks[$i]}"
      input_folder="${input_folders[$i]}"
      mkdir -p ./${bb}_regenie
      for pheno in "${phenotypes[@]}"; do
        gcloud storage cp ${input_folder}/*_chr~{chr}_*_${pheno}.regenie.gz ./${bb}_regenie/
      done
    done

    # ── Resolve annotation source ────────────────────────────────────
    # Exactly one of schema_directory or annotation_manifest should be non-empty, if aggs run.
    schema_dir_input="~{schema_directory}"
    schema_dir_input="${schema_dir_input//\{chr\}/~{chr}}"
    annotation_manifest="~{annotation_manifest}"
    use_parquet=false
    acat_annot_flags=""
    ann_args=""

    if [[ -n "${schema_dir_input}" && -n "${annotation_manifest}" ]]; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: both schema_directory and annotation_manifest provided, they cannot be combined" >> ${log_file}
      exit 1
    fi

    if [[ -n "${schema_dir_input}" ]]; then
      # Parquet schema directory — download all 6 parquet files
      use_parquet=true
      mkdir -p annots/schema
      gcloud storage cp "${schema_dir_input}/*.parquet" ./annots/schema/
      acat_annot_flags="--schema-dir ./annots/schema"
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] Using parquet schema from ${schema_dir_input}" >> ${log_file}

    elif [[ -n "${annotation_manifest}" ]]; then
      # Manifest string: "name,mask,annot,set;name2,mask2,annot2,set2;..."
      # Parse into parallel arrays and download mask + annotation files per chr.
      mkdir -p annots
      masks_args=()
      annots_args=()

      IFS=';' read -ra triplets <<< "${annotation_manifest}"
      for triplet in "${triplets[@]}"; do
        IFS=',' read -r entry_name mask_file annot_file _setlist <<< "${triplet}"
        # Download mask file (chr-independent, or with {chr} placeholder)
        mask_resolved="${mask_file//\{chr\}/~{chr}}"
        gcloud storage cp "${mask_resolved}" ./annots/
        local_mask="./annots/$(basename "${mask_resolved}")"
        # Decompress if needed
        if [[ "${local_mask}" == *.zst ]]; then
          zstd -d --rm "${local_mask}"
          local_mask="${local_mask%.zst}"
        fi
        masks_args+=("${local_mask}")

        # Download annotation file (typically has {chr} placeholder)
        annot_resolved="${annot_file//\{chr\}/~{chr}}"
        gcloud storage cp "${annot_resolved}" ./annots/
        local_annot="./annots/$(basename "${annot_resolved}")"
        if [[ "${local_annot}" == *.zst ]]; then
          zstd -d --rm "${local_annot}"
          local_annot="${local_annot%.zst}"
        fi
        annots_args+=("${local_annot}")
      done

      acat_annot_flags="--masks-tsvs ${masks_args[*]} --annotations-tsvs ${annots_args[*]}"
      ann_args="${annots_args[*]}"
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] Using manifest annotations: ${#triplets[@]} entries" >> ${log_file}

    elif [[ "${do_acat}" == "true" || "${do_burden}" == "true" ]]; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: neither schema_directory nor annotation_manifest provided" >> ${log_file}
      exit 1
    fi

    # ── Build cohort file for fedcat meta ────────────────────────────
    # Format: "# cohort_name" header lines followed by file paths.
    cohort_file="cohort_files_chr~{chr}.txt"
    for pair in ~{sep=' ' stratum_biobank_pairs}; do
      stratum="${pair%%_*}"
      bb="${pair##*_}"
      echo "# ${pair}" >> "${cohort_file}"
      for f in ./${bb}_regenie/*_${stratum}_chr~{chr}_*.regenie.gz; do
        [[ -e "$f" ]] && echo "$f" >> "${cohort_file}"
      done
    done

    # ── Meta-analysis per result type ─────────────────────────────────
    mkdir -p meta_outputs
    ext=".zst"

    run_meta_for_type() {
      local res_type="$1"
      local should_run="$2"
      [[ "${should_run}" != "true" ]] && return 0

      for pheno in "${phenotypes[@]}"; do
        # Build a per-phenotype cohort file filtered to this res_type
        local pheno_cohort="cohort_${pheno}_${res_type}_chr~{chr}.txt"
        for pair in ~{sep=' ' stratum_biobank_pairs}; do
          stratum="${pair%%_*}"
          bb="${pair##*_}"
          echo "# ${pair}" >> "${pheno_cohort}"
          if [[ "${res_type}" == "single_variants" ]]; then
            pattern="./${bb}_regenie/*_${stratum}_chr~{chr}*_singlevar_${pheno}.regenie.gz"
          elif [[ "${res_type}" == "single_variants_unconditional" ]]; then
            pattern="./${bb}_regenie/*_${stratum}_chr~{chr}*_singlevar_unconditional_${pheno}.regenie.gz"
          elif [[ "${res_type}" == "burden" ]]; then
            pattern="./${bb}_regenie/*_${stratum}_chr~{chr}_*_burden_${pheno}.regenie.gz"
          elif [[ "${res_type}" == "mac10" ]]; then
            pattern="./${bb}_regenie/*_${stratum}_chr~{chr}_*_mac10_${pheno}.regenie.gz"
          fi
          for f in ${pattern}; do
            [[ -e "$f" ]] && echo "$f" >> "${pheno_cohort}"
          done
        done

        local out_path="meta_outputs/${pheno}_${res_type}_meta_chr~{chr}.tsv${ext}"
        fedcat meta \
          --cohort-file "${pheno_cohort}" \
          --output-file "${out_path}" \
          --res-type "${res_type}" \
          --compression zst \
          2>&1 | tee -a ${log_file}
      done
    }

    run_meta_for_type "single_variants" "${do_cond_singles}"
    run_meta_for_type "single_variants_unconditional" "${do_uncond_singles}"
    run_meta_for_type "mac10" "${do_mac10}"
    run_meta_for_type "burden" "${do_burden}"

    # ACAT singles: a stratum enters each AAF bin only where the variant is rare in
    # it. Its presence in the output folder is what marks the ACAT as built this way.
    if [[ "${do_acat}" == "true" && "~{acat_singles_stratum_rare}" == "true" ]]; then
      if [[ "${aggs_uncond}" == "true" ]]; then singles_type="single_variants_unconditional"; else singles_type="single_variants"; fi
      for pheno in "${phenotypes[@]}"; do
        singles_cohort="cohort_${pheno}_${singles_type}_chr~{chr}.txt"
        [[ -f "${singles_cohort}" ]] || continue
        fedcat meta \
          --cohort-file "${singles_cohort}" \
          --output-file "meta_outputs/${pheno}_acat_singles_meta_chr~{chr}.tsv${ext}" \
          --res-type "${singles_type}" \
          --aaf-bins "~{aaf_bins}" \
          --compression zst \
          2>&1 | tee -a ${log_file}
      done
    fi

    # Upload meta-analysis outputs (the log is uploaded last, at task end)
    shopt -s nullglob
    meta_zst=(./meta_outputs/*.tsv.zst)
    if (( ${#meta_zst[@]} )); then
      gcloud storage cp "${meta_zst[@]}" ~{output_folder}/
    fi

    # ── ACAT ──────────────────────────────────────────────────────────
    if [[ "${do_acat}" == "true" ]]; then
      for pheno in "${phenotypes[@]}"; do
        if [[ "~{acat_singles_stratum_rare}" == "true" ]]; then
          singles="meta_outputs/${pheno}_acat_singles_meta_chr~{chr}.tsv${ext}"
        elif [[ "${aggs_uncond}" == "true" ]]; then
          singles="meta_outputs/${pheno}_single_variants_unconditional_meta_chr~{chr}.tsv${ext}"
        else
          singles="meta_outputs/${pheno}_single_variants_meta_chr~{chr}.tsv${ext}"
        fi
        mac10_file="meta_outputs/${pheno}_mac10_meta_chr~{chr}.tsv${ext}"
        acat_out="${pheno}_acat_chr~{chr}.tsv${ext}"

        if [[ -f "${singles}" ]]; then
          # mac10 burden is optional; without it ACAT combines the single variants alone
          acat_burden_flag=""
          if [[ -f "${mac10_file}" ]]; then acat_burden_flag="--burden-file ${mac10_file}"; fi
          fedcat acat \
            --single-variants-file "${singles}" \
            ${acat_burden_flag} \
            ${acat_annot_flags} \
            --output-file "${acat_out}" \
            --n-threads ~{n_cpu - 1} \
            --aggregate-freqs "~{aaf_bins}" \
            ~{if acat_full_bounds then "" else "--no-full-bounds"} \
            2>&1 | tee -a ${log_file}
        else
          echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skipping ACAT for ${pheno} chr~{chr}: missing meta outputs" >> ${log_file}
        fi
      done
    fi

    # Upload ACAT outputs (the log is uploaded last, at task end)
    shopt -s nullglob
    acat_zst=(./*_acat_chr~{chr}.tsv.zst)
    if (( ${#acat_zst[@]} )); then
      gcloud storage cp "${acat_zst[@]}" ~{output_folder}/
    fi

    # Step-0 subsetting
    if [[ ~{create_step0_subset} == true ]]; then
      # Pass annotation source to step0 subsetting
      step0_ann_flags=""
      if [[ "${use_parquet}" == "true" ]]; then
        step0_ann_flags="--schema-dir ./annots/schema"
      elif [[ -n "${ann_args}" ]]; then
        step0_ann_flags="--annotations ${ann_args} --masks-tsvs ${masks_args[*]}"
      fi

      for pheno in "${phenotypes[@]}"; do
        step0_file_flags=""
        if [[ "${do_acat}" == "true" ]]; then
          acat_f="${pheno}_acat_chr~{chr}.tsv.zst"
          [[ -f "${acat_f}" ]] && step0_file_flags+=" --acat-file ${acat_f}"
        fi
        if [[ "${do_burden}" == "true" ]]; then
          burden_f="meta_outputs/${pheno}_burden_meta_chr~{chr}.tsv.zst"
          [[ -f "${burden_f}" ]] && step0_file_flags+=" --burden-file ${burden_f}"
        fi
        if [[ "${do_uncond_singles}" == "true" ]]; then
          singles_f="meta_outputs/${pheno}_single_variants_unconditional_meta_chr~{chr}.tsv.zst"
          [[ -f "${singles_f}" ]] && step0_file_flags+=" --singles-file ${singles_f}"
        fi

        if [[ -n "${step0_file_flags}" ]]; then
          fedcat subset-results \
            --phenotypes "${pheno}" \
            --chr ~{chr} \
            --singles-p-threshold ~{cond_p_singles_retest} \
            --agg-p-threshold ~{cond_p_agg_retest} \
            ${step0_ann_flags} \
            ${step0_file_flags} \
            2>&1 | tee -a ${log_file}
        fi
      done

      # Each type goes to its own folder, with the marker the step0_subset phase checks for
      upload_step0() {
        local cond_type="$1" folder="$2"
        local tars=(./*_${cond_type}_cond_step0_subset_chr~{chr}.tar.gz)
        if (( ${#tars[@]} )); then
          gcloud storage cp "${tars[@]}" "${folder}/"
        fi
        local marker="~{run_prefix}_step0_subset_chr~{chr}.log"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${cond_type} step0 subset built by the meta-analysis task" > "${marker}"
        gcloud storage cp "${marker}" "${folder}/"
      }
      shopt -s nullglob
      singles_folder="~{step0_singles_folder}"
      agg_folder="~{step0_agg_folder}"
      if [[ "${do_uncond_singles}" == "true" ]]; then
        upload_step0 singles "${singles_folder:-~{output_folder}}"
      fi
      if [[ "${do_burden}" == "true" || "${do_acat}" == "true" ]]; then
        upload_step0 agg "${agg_folder:-~{output_folder}}"
      fi
    fi

    # upload the log last so its presence marks the run as complete (see the check at the top)
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Meta-analysis task completed" >> ${log_file}
    gcloud storage cp "${log_file}" ~{output_folder}/
  >>>

  runtime {
    docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/hwright0/fedcat_meta:0.8.1"
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + disk
    maxRetries: max_retries
    zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
    cpuPlatform: "Intel Cascade Lake"
    preemptible: preemptible_tries
  }
}
