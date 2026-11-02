version 1.0

workflow regenie_step2_aggregates_aou {
    input {
        File phenotype_file
        File covariates_file
        String conditional_snps_path = ""
        String trait_Binary_or_Continuous
        String pheno_transformation_RINT_or_MCC_or_RAW
        String comma_seperated_list_of_categorical_covariates_or_NONE
        # One entry per phenotype group: its comma-separated phenotypes, the step1
        # prefixes holding their predictions (comma-separated), a tag for its output
        # names ("" for none) and its tasks as [stratum, chr, start, end, exclcovars]
        Array[String] pheno_groups
        Array[String] step1_run_prefixes
        Array[String] group_labels
        Array[Array[Array[String]]] units
        String bucket
        String run_prefix

        String stratify_covar = ""
        String stratify_value = ""
        Int n_cpu = 8
        Int mem_gb = 64
        String disk = "1000 HDD"
        Int max_retries = 0
        Int preemptible_tries = 3
        String is_step1_done = "yes"
        String step1_directory = "regenie_step1_output"
        String step1_input_template
        String results_directory = "wgs_analysis2/v2026.2"
        String pfile_template = "wgs_v8_pgens_genoqc/chr{chr}_genoqc_trimmed_merged"
        String annotation_manifest
        Boolean coding_only = false
        String aaf_bins = "0.001"
        String regenie_version = "4.1.3"
        String dependencies_path = "gs://fc-secure-ab157c8b-aa91-4d13-b17f-8276c7b7fbcb/workflows/regenie_wgs_workflow/v1.7.1/dependencies"
    }

    scatter(idx in range(length(pheno_groups))){
        scatter(unit in units[idx]){
            call regenie_step2_aggregates {
                input:
                    phenotype_tsv = phenotype_file,
                    covariate_tsv = covariates_file,
                    exclcovars = unit[4],
                    conditional = conditional_snps_path,
                    type = trait_Binary_or_Continuous,
                    invnorm = pheno_transformation_RINT_or_MCC_or_RAW,
                    catvars = comma_seperated_list_of_categorical_covariates_or_NONE,
                    chr = unit[1],
                    chrom_start = unit[2],
                    chrom_end = unit[3],
                    stratum = unit[0],
                    phenos_to_run = pheno_groups[idx],
                    bucket = bucket,
                    run_prefix = run_prefix,
                    step1_run_prefix = step1_run_prefixes[idx],
                    pheno_group = group_labels[idx],
                    stratify_covar = stratify_covar,
                    stratify_value = stratify_value,
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
                    coding_only = coding_only,
                    aaf_bins = aaf_bins,
                    regenie_version = regenie_version,
                    dependencies_path = dependencies_path
            }
        }
    }
    output {
        Array[Array[String]] step2_aggregate_completions = regenie_step2_aggregates.step2_aggregates_completed
        String step2_aggregates_done = "YES"
    }
}

task regenie_step2_aggregates {

    input {
        File phenotype_tsv
        File covariate_tsv
        String exclcovars
        String conditional
        String type
        String invnorm
        String catvars
        String chr
        String chrom_start
        String chrom_end
        String stratum
        String phenos_to_run
        String bucket
        String run_prefix
        String step1_run_prefix
        String pheno_group
        String stratify_covar = ""
        String stratify_value = ""
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
        Boolean coding_only
        String aaf_bins
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

        ./regenie_step2_aggregates_aou.sh \
            "~{phenotype_tsv}" \
            "~{covariate_tsv}" \
            "~{exclcovars}" \
            "~{conditional}" \
            "~{type}" \
            "~{invnorm}" \
            "~{catvars}" \
            "~{chr}" \
            "~{chrom_start}" \
            "~{chrom_end}" \
            "~{stratum}" \
            "~{phenos_to_run}" \
            "~{bucket}" \
            "~{run_prefix}" \
            "~{step1_directory}" \
            "~{step1_input_template}" \
            "~{results_directory}" \
            "~{pfile_template}" \
            "~{annotation_manifest}" \
            "~{coding_only}" \
            "~{aaf_bins}" \
            "~{step1_run_prefix}" \
            "~{stratify_covar}" \
            "~{stratify_value}" \
            "~{pheno_group}"

    >>>

    output {
        String step2_aggregates_completed = "chr~{chr}_~{chrom_start}_~{chrom_end}_~{stratum}_burden_completed"
    }

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
