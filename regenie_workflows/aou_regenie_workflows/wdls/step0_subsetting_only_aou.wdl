version 1.0

workflow step0_subsetting_only_aou {

  input {
    Array[String] chromosomes
    # meta_folder holds the meta-analysis outputs; output_folder receives this threshold's subsets
    String meta_folder
    String output_folder
    String run_prefix
    Array[String] phenotypes
    # Annotation source — provide exactly ONE of these two:
    #   schema_directory  → full GCS path to pre-built optimised parquet schema
    #   annotation_manifest → semicolon-separated AnnotationEntry triplets:
    #       name,mask_file,annot_file,setlist_file;name2,...
    #       Paths should be full GCS paths with {chr} placeholders.
    String schema_directory = ""
    String annotation_manifest = ""
    String dependencies_path = "gs://fc-secure-ab157c8b-aa91-4d13-b17f-8276c7b7fbcb/workflows/regenie_wgs_workflow/v1.7.1/dependencies"

    # Which meta outputs to subset; they must exist in meta_folder.
    Boolean run_uncond_singles = true
    Boolean run_burden = true
    Boolean run_acat = true
    Boolean acat_singles_stratum_rare = true
    
    Float cond_p_singles_retest = 0.000001
    Float cond_p_agg_retest = 0.000001

    Int n_cpu = 8
    Int mem_gb = 32
    String disk = "400 HDD"
    Int max_retries = 0
    Int preemptible_tries = 0
  }

  scatter (chr in chromosomes) {
    call subset_per_chromosome {
      input:
        chr = chr,
        meta_folder = meta_folder,
        output_folder = output_folder,
        run_prefix = run_prefix,
        phenotypes = phenotypes,
        schema_directory = schema_directory,
        annotation_manifest = annotation_manifest,
        dependencies_path = dependencies_path,
        run_uncond_singles = run_uncond_singles,
        run_burden = run_burden,
        run_acat = run_acat,
        acat_singles_stratum_rare = acat_singles_stratum_rare,
        cond_p_singles_retest = cond_p_singles_retest,
        cond_p_agg_retest = cond_p_agg_retest,
        n_cpu = n_cpu,
        mem_gb = mem_gb,
        disk = disk,
        max_retries = max_retries,
        preemptible_tries = preemptible_tries
    }
  }
}

task subset_per_chromosome {
  input {
    String chr
    String meta_folder
    String output_folder
    String run_prefix
    Array[String] phenotypes
    String schema_directory
    String annotation_manifest
    String dependencies_path
    Boolean run_uncond_singles
    Boolean run_burden
    Boolean run_acat
    Boolean acat_singles_stratum_rare
    Float cond_p_singles_retest
    Float cond_p_agg_retest
    Int n_cpu
    Int mem_gb
    String disk
    Int max_retries
    Int preemptible_tries
  }

  command <<<
    set -euox pipefail
    version="1.7.1"
    log_file=~{run_prefix}_step0_subset_chr~{chr}.log

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting step0 subsetting-only task, version ${version}" > ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] chr: ~{chr}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] meta_folder: ~{meta_folder}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] output_folder: ~{output_folder}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] phenotypes: ~{sep="," phenotypes}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] schema_directory: ~{schema_directory}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] annotation_manifest: ~{annotation_manifest}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] run_uncond_singles: ~{run_uncond_singles}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] run_burden: ~{run_burden}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] run_acat: ~{run_acat}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] cond_p_singles_retest: ~{cond_p_singles_retest}" >> ${log_file}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] cond_p_agg_retest: ~{cond_p_agg_retest}" >> ${log_file}

    # Skip if this chromosome already completed. We upload the log last, so its
    # presence in the output folder means the run finished (independent of caching).
    if gcloud storage objects describe "~{output_folder}/${log_file}" >/dev/null 2>&1; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] Log file ${log_file} already exists, assuming complete and exiting"
      exit 0
    fi

    # Pull in scripts and install fedcat from the dependency path
    gcloud config set storage/parallel_composite_upload_enabled False
    gcloud storage cp ~{dependencies_path}/* .
    chmod 777 *

    # ── Resolve annotation source ────────────────────────────────────
    schema_dir_input="~{schema_directory}"
    schema_dir_input="${schema_dir_input//\{chr\}/~{chr}}"
    annotation_manifest="~{annotation_manifest}"
    use_parquet=false
    ann_args=""

    if [[ -n "${schema_dir_input}" && -n "${annotation_manifest}" ]]; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: both schema_directory and annotation_manifest provided, they cannot be combined" >> ${log_file}
      exit 1
    fi

    if [[ -n "${schema_dir_input}" ]]; then
      use_parquet=true
      mkdir -p annots/schema
      gcloud storage cp "${schema_dir_input}/*.parquet" ./annots/schema/
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] Using parquet schema from ${schema_dir_input}" >> ${log_file}

    elif [[ -n "${annotation_manifest}" ]]; then
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

        annot_resolved="${annot_file//\{chr\}/~{chr}}"
        gcloud storage cp "${annot_resolved}" ./annots/
        local_annot="./annots/$(basename "${annot_resolved}")"
        if [[ "${local_annot}" == *.zst ]]; then
          zstd -d --rm "${local_annot}"
          local_annot="${local_annot%.zst}"
        fi
        annots_args+=("${local_annot}")
      done

      ann_args="${annots_args[*]}"
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] Using manifest annotations: ${#triplets[@]} entries" >> ${log_file}

    elif [[ "~{run_acat}" == "true" || "~{run_burden}" == "true" ]]; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: neither schema_directory nor annotation_manifest provided" >> ${log_file}
      exit 1
    fi

    # Build step0 annotation flags
    step0_ann_flags=""
    if [[ "${use_parquet}" == "true" ]]; then
      step0_ann_flags="--schema-dir ./annots/schema"
    elif [[ -n "${ann_args}" ]]; then
      step0_ann_flags="--annotations ${ann_args} --masks-tsvs ${masks_args[*]}"
    fi

    # Download fedcat meta-analysis outputs (already produced) for this chromosome.
    # They are read as .tsv.zst directly.
    phenotypes=(~{sep=' ' phenotypes})
    mkdir -p meta_outputs

    for phenotype in "${phenotypes[@]}"; do
      if [[ "~{run_uncond_singles}" == "true" ]]; then
        gcloud storage cp "~{meta_folder}/${phenotype}_single_variants_unconditional_meta_chr~{chr}.tsv.zst" ./meta_outputs/
      fi
      if [[ "~{run_burden}" == "true" ]]; then
        gcloud storage cp "~{meta_folder}/${phenotype}_burden_meta_chr~{chr}.tsv.zst" ./meta_outputs/
      fi
      if [[ "~{run_acat}" == "true" ]]; then
        gcloud storage cp "~{meta_folder}/${phenotype}_acat_chr~{chr}.tsv.zst" ./
        # Conditioning recomputes ACAT each step with this rule, so step 0 must match
        acat_singles_found=false
        if gcloud storage ls "~{meta_folder}/${phenotype}_acat_singles_meta_chr~{chr}.tsv.zst" > /dev/null 2>&1; then
          acat_singles_found=true
        fi
        if [[ "${acat_singles_found}" != "~{acat_singles_stratum_rare}" ]]; then
          echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: the meta-analysis ACAT for ${phenotype} chr~{chr} was built with acat_singles_stratum_rare=${acat_singles_found}, but this run uses ~{acat_singles_stratum_rare}; re-run the meta-analysis or change the config" >> ${log_file}
          exit 1
        fi
      fi
    done

    # Run the subsetting per phenotype; tarballs written at top-level.
    for phenotype in "${phenotypes[@]}"; do
      step0_file_flags=""
      if [[ "~{run_acat}" == "true" ]]; then
        acat_f="${phenotype}_acat_chr~{chr}.tsv.zst"
        [[ -f "${acat_f}" ]] && step0_file_flags+=" --acat-file ${acat_f}"
      fi
      if [[ "~{run_burden}" == "true" ]]; then
        burden_f="meta_outputs/${phenotype}_burden_meta_chr~{chr}.tsv.zst"
        [[ -f "${burden_f}" ]] && step0_file_flags+=" --burden-file ${burden_f}"
      fi
      if [[ "~{run_uncond_singles}" == "true" ]]; then
        singles_f="meta_outputs/${phenotype}_single_variants_unconditional_meta_chr~{chr}.tsv.zst"
        [[ -f "${singles_f}" ]] && step0_file_flags+=" --singles-file ${singles_f}"
      fi

      if [[ -n "${step0_file_flags}" ]]; then
        fedcat subset-results \
          --phenotypes "${phenotype}" \
          --chr ~{chr} \
          --singles-p-threshold ~{cond_p_singles_retest} \
          --agg-p-threshold ~{cond_p_agg_retest} \
          ${step0_ann_flags} \
          ${step0_file_flags} \
          2>&1 | tee -a ${log_file}
      fi
    done

    shopt -s nullglob
    subset_tars=(./*_step0_subset_chr~{chr}.tar.gz)
    if (( ${#subset_tars[@]} )); then
      gcloud storage cp "${subset_tars[@]}" ~{output_folder}/
    fi

    # upload the log last so its presence marks the run as complete (see the check at the top)
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step0 subsetting task completed" >> ${log_file}
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
