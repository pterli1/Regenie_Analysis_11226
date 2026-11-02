version 1.0

workflow regenie_step1_scatter_aou {
    input {
        File phenotype_file
        File covariates_file
        String trait_Binary_or_Continuous
        String pheno_transformation_RINT_or_MCC_or_RAW
        String comma_seperated_list_of_categorical_covariates_or_NONE
        # regenie --covarExcludeList for the fit, or "None"
        String exclcovars
        # One entry per phenotype set: its comma-separated phenotypes, output
        # prefix and the strata still to fit
        Array[String] pheno_sets
        Array[String] run_prefixes
        Array[Array[String]] strata_by_set
        String bucket
        String stratify_covar = ""
        String stratify_value = ""
        Int n_cpu = 8
        Int mem_gb = 32
        String disk = "350 SSD"
        Int max_retries = 0
        Int preemptible_tries = 1
        String step1_directory = "regenie_step1_output"
        String step1_input_template
        String regenie_version = "4.1.3"
        String dependencies_path = "gs://fc-secure-ab157c8b-aa91-4d13-b17f-8276c7b7fbcb/workflows/regenie_wgs_workflow/v1.7.1/dependencies"
    }

    scatter(set_idx in range(length(pheno_sets))){
        scatter(stratum in strata_by_set[set_idx]){
            call run_regenie_step1 {
            input:
                phenotype_tsv = phenotype_file,
                covariate_tsv = covariates_file,
                type = trait_Binary_or_Continuous,
                invnorm = pheno_transformation_RINT_or_MCC_or_RAW,
                catvars = comma_seperated_list_of_categorical_covariates_or_NONE,
                exclcovars = exclcovars,
                stratum = stratum,
                phenos_to_run = pheno_sets[set_idx],
                bucket = bucket,
                run_prefix = run_prefixes[set_idx],
                stratify_covar = stratify_covar,
                stratify_value = stratify_value,
                n_cpu = n_cpu,
                mem_gb = mem_gb,
                disk = disk,
                max_retries = max_retries,
                preemptible_tries = preemptible_tries,
                step1_directory = step1_directory,
                step1_input_template = step1_input_template,
                regenie_version = regenie_version,
                dependencies_path = dependencies_path
            }
        }
    }

    output{
        String step1_done = "yes"
    }
}

task run_regenie_step1 {

    input {
        File phenotype_tsv
        File covariate_tsv
        String type
        String invnorm
        String catvars
        String exclcovars
        String stratum
        String phenos_to_run
        String bucket
        String run_prefix
        String stratify_covar = ""
        String stratify_value = ""
        Int n_cpu
        Int mem_gb
        String disk
        Int max_retries
        Int preemptible_tries
        String step1_directory
        String step1_input_template
        String regenie_version
        String dependencies_path
    }

    command <<<
        set -euxo pipefail

        version="1.7.1"
        log_file=~{run_prefix}_~{type}_~{invnorm}_Step1_~{stratum}.workflow.log
        bucket=~{bucket}
        step1_directory=~{step1_directory}
        step1_input_template=~{step1_input_template}

        # the log is uploaded last, so its presence means this fit completed
        if gcloud storage objects describe "${bucket}/${step1_directory}/${log_file}" >/dev/null 2>&1; then
            echo "Log file ${log_file} already exists, assuming complete and exiting"
            exit 0
        fi

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting regenie step1 task, version ${version}" > ${log_file}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] phenotype_tsv: ~{phenotype_tsv}" >> ${log_file}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] covariate_tsv: ~{covariate_tsv}" >> ${log_file}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] type: ~{type}" >> ${log_file}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] invnorm: ~{invnorm}" >> ${log_file}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] catvars: ~{catvars}" >> ${log_file}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] exclcovars: ~{exclcovars}" >> ${log_file}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] stratum: ~{stratum}" >> ${log_file}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] phenos_to_run: ~{phenos_to_run}" >> ${log_file}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] run_prefix: ~{run_prefix}" >> ${log_file}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] step1_directory: ~{step1_directory}" >> ${log_file}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] stratify_covar: ~{stratify_covar}" >> ${log_file}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] stratify_value: ~{stratify_value}" >> ${log_file}

        # copy in scripts
        gcloud storage cp ~{dependencies_path}/* .
        # set permissions
        chmod 777 *
        # point regenie at the requested build; --version fails here if the image lacks it
        ln -sf /usr/local/bin/regenie-~{regenie_version} /usr/local/bin/regenie
        regenie --version
        source shared_functions.sh

        # resolve step1 template into step1_input_path/_directory/_name
        split_template "${step1_input_template}" stratum "~{stratum}" step1_input

        # decompress phenotype/covariate files (handles .gz/.zst/plain)
        decompress_to "~{covariate_tsv}" covariate_tsv
        decompress_to "~{phenotype_tsv}" phenotype_tsv

        # download step 1 input files
        gcloud storage cp "${step1_input_directory}/${step1_input_name}*" .
        # regenie reads this fileset directly below, so normalise the .fam first
        strip_fam_header "${step1_input_name}"

        # subset covariate and phenotype files by stratum
        extract_stratum_ids "${step1_input_name}" > stratum_ids
        subset_to_ids stratum_ids covariate_tsv
        subset_to_ids stratum_ids phenotype_tsv
        # then by covariate value, if this is a covariate-stratified run
        subset_to_covar_value "~{stratify_covar}" "~{stratify_value}" covariate_tsv phenotype_tsv
        build_subset_ids covariate_tsv phenotype_tsv > subset_ids

        # filter step1 input to sample mac >= 10
        plink_com="plink2 \
            --keep subset_ids \
            --mac 10 \
            --make-pgen \
            --out ${step1_input_name}_mac10"

        # pgen vs bed input
        if [[ -f ${step1_input_name}.pgen ]]; then
            plink_com=${plink_com}" --pfile ${step1_input_name}"
        else
            plink_com=${plink_com}" --bfile ${step1_input_name}"
        fi

        $plink_com

        # build command
        regenie_com="regenie \
            --step 1 \
            --pgen ${step1_input_name}_mac10 \
            --phenoFile phenotype_tsv \
            --covarFile covariate_tsv \
            --phenoColList ~{phenos_to_run} \
            --bsize 1000 \
            --gz \
            --lowmem \
            --threads ~{n_cpu} \
            --out ~{run_prefix}_~{type}_~{invnorm}_Step1_~{stratum}"
        
        # add other args as required
        # categorical covariates
        if [[ ~{catvars} != "NONE" && ~{catvars} != "" ]]; then
            regenie_com=${regenie_com}" \
            --catCovarList ~{catvars} \
            --maxCatLevels 30"
        fi
        # covariates excluded from every chromosome
        if [[ "~{exclcovars}" != "None" ]]; then
            regenie_com=${regenie_com}" --covarExcludeList ~{exclcovars}"
        fi
        # trait type/transform options
        if [ ~{type} == "Continuous" ]; then
            if [ "~{invnorm}" == "RINT" ]; then
                regenie_com=${regenie_com}" --apply-rint"
            elif [ "~{invnorm}" == "MCC" ]; then
                regenie_com=${regenie_com}" --mcc  --mcc-thr 0.01"
            fi
        else
            regenie_com=${regenie_com}" --bt"
        fi

        # run regenie
        $regenie_com

        # upload output
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed regenie step1 task" >> ${log_file}
        mkdir step1_out
        mv ~{run_prefix}_~{type}_~{invnorm}_Step1_~{stratum}* step1_out/
        mv step1_out/${log_file} .
        gcloud storage cp step1_out/* ${bucket}/${step1_directory}/
        gcloud storage cp ${log_file} ${bucket}/${step1_directory}/
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
