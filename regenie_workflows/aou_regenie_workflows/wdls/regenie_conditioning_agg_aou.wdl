version 1.0

workflow regenie_conditioning_agg_aou {
    input {
        File phenotype_file
        File covariates_file
        String conditional_snps_path = ""
        String trait_Binary_or_Continuous
        String pheno_transformation_RINT_or_MCC_or_RAW
        String comma_seperated_list_of_categorical_covariates_or_NONE
        Array[String] strata_array
        Array[Array[Array[String]]] conditioning_blocks
        String bucket
        String run_prefix
        String step1_run_prefix = run_prefix
        String stratify_covar = ""
        String stratify_value = ""
        String run_suffix
        String pheno_to_run
        Int n_cpu = 8
        Int mem_gb = 64
        String disk = "1000 HDD"
        Int max_retries = 0
        Int preemptible_tries = 3
        String step1_directory = "regenie_step1_output"
        String step1_input_template
        String results_directory = "wgs_analysis2/v2026.2"
        String pfile_template = "wgs_v8_pgens_genoqc/chr{chr}_genoqc_trimmed_merged"
        String annotation_manifest
        String aaf_bins = "0.001"
        Int min_arm_allele_count = 0
        String cond_engine = "regenie"
        Float cond_rsq_threshold = 0.9
        Boolean acat_singles_stratum_rare = false
        String cond_mancojo_effect_mode = "imputeNA"
        Int cond_step
        String regenie_version = "4.1.3"
        String dependencies_path = "gs://fc-secure-ab157c8b-aa91-4d13-b17f-8276c7b7fbcb/workflows/regenie_wgs_workflow/v1.7.1/dependencies"
    }

    scatter(idx in range(length(strata_array))){
        Array[Array[String]] blocks = conditioning_blocks[idx]
        scatter(block in blocks){
            call regenie_conditioning_agg {
                input:
                    phenotype_tsv = phenotype_file,
                    covariate_tsv = covariates_file,
                    exclcovars = block[1],
                    conditional = conditional_snps_path,
                    type = trait_Binary_or_Continuous,
                    invnorm = pheno_transformation_RINT_or_MCC_or_RAW,
                    catvars = comma_seperated_list_of_categorical_covariates_or_NONE,
                    chr = block[0],
                    aggregates = block[3],
                    cond_aggregates = block[2],
                    stratum = strata_array[idx],
                    bucket = bucket,
                    run_prefix = run_prefix,
                    step1_run_prefix = step1_run_prefix,
                    stratify_covar = stratify_covar,
                    stratify_value = stratify_value,
                    run_suffix = run_suffix,
                    pheno_to_run = pheno_to_run,
                    n_cpu = n_cpu,
                    mem_gb = mem_gb,
                    disk = disk,
                    max_retries = max_retries,
                    preemptible_tries = preemptible_tries,
                    step1_directory = step1_directory,
                    step1_input_template = step1_input_template,
                    results_directory = results_directory,
                    pfile_template = pfile_template,
                    annotation_manifest = annotation_manifest,
                    aaf_bins = aaf_bins,
                    min_arm_allele_count = min_arm_allele_count,
                    cond_engine = cond_engine,
                    cond_rsq_threshold = cond_rsq_threshold,
                    acat_singles_stratum_rare = acat_singles_stratum_rare,
                    cond_mancojo_effect_mode = cond_mancojo_effect_mode,
                    cond_step = cond_step,
                    regenie_version = regenie_version,
                    dependencies_path = dependencies_path
            }
        }
    }
}

task regenie_conditioning_agg {

    input {
        File phenotype_tsv
        File covariate_tsv
        String exclcovars
        String conditional
        String type
        String invnorm
        String catvars
        String chr
        String aggregates
        String cond_aggregates
        String stratum
        String bucket
        String run_prefix
        String step1_run_prefix
        String stratify_covar = ""
        String stratify_value = ""
        String run_suffix
        String pheno_to_run
        Int n_cpu
        Int mem_gb
        String disk
        Int max_retries
        Int preemptible_tries
        String step1_directory
        String step1_input_template
        String results_directory
        String pfile_template
        String annotation_manifest
        String aaf_bins
        Int min_arm_allele_count
        String cond_engine
        Float cond_rsq_threshold
        Boolean acat_singles_stratum_rare
        String cond_mancojo_effect_mode
        Int cond_step
        String regenie_version
        String dependencies_path
    }

    command <<<
        set -euxo pipefail
        # copy in scripts and binaries
        gcloud storage cp ~{dependencies_path}/* .
        # set permissions
        chmod 777 *
        # point regenie at the requested build; --version fails here if the image lacks it
        ln -sf /usr/local/bin/regenie-~{regenie_version} /usr/local/bin/regenie
        regenie --version

        ./regenie_conditioning_agg_aou.sh \
            "~{phenotype_tsv}" \
            "~{covariate_tsv}" \
            "~{exclcovars}" \
            "~{conditional}" \
            "~{type}" \
            "~{invnorm}" \
            "~{catvars}" \
            "~{chr}" \
            "~{aggregates}" \
            "~{cond_aggregates}" \
            "~{stratum}" \
            "~{bucket}" \
            "~{run_prefix}" \
            "~{run_suffix}" \
            "~{pheno_to_run}" \
            "~{step1_directory}" \
            "~{step1_input_template}" \
            "~{results_directory}" \
            "~{pfile_template}" \
            "~{annotation_manifest}" \
            "~{aaf_bins}" \
            "~{step1_run_prefix}" \
            "~{stratify_covar}" \
            "~{stratify_value}" \
            "~{min_arm_allele_count}" \
            "~{cond_engine}" \
            "~{cond_rsq_threshold}" \
            "~{acat_singles_stratum_rare}" \
            "~{cond_mancojo_effect_mode}" \
            "~{cond_step}"
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
