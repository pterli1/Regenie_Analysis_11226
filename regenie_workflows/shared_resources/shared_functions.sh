#!/bin/bash
# Shared helper functions for the regenie step scripts (UKB applets and AoU
# workflows). Source this file, then call the helpers; they are not meant to be
# run directly.

# ── Conditional decompression ────────────────────────────────────────────────
# Materialise $1 (a local file that may be gzip/zstd compressed or plain) as the
# uncompressed file $2. The compression is detected from the source extension so
# callers no longer have to assume a format.
decompress_to() {
  local src=$1 dest=$2
  case "$src" in
    *.gz)
      gunzip -c "$src" > "$dest"
      ;;
    *.zst)
      # zstd CLI if present, otherwise the bundled plink2 (always on PATH here).
      if command -v zstd >/dev/null 2>&1; then
        zstd -dc "$src" > "$dest"
      else
        plink2 --zst-decompress "$src" "$dest"
      fi
      ;;
    *)
      [[ "$src" == "$dest" ]] || cp "$src" "$dest"
      ;;
  esac
}

# ── Template resolution ──────────────────────────────────────────────────────
# Substitute {$2} -> $3 in template $1, then split the result into
# <$4>_path / <$4>_directory / <$4>_name in the caller's scope.
#   split_template "$pfile_template" chr "$chr" pfile   -> pfile_path/_directory/_name
#   split_template "$step1_input_template" stratum "$stratum" step1_input
split_template() {
  local template=$1 placeholder=$2 value=$3 out=$4 path
  path=$(echo "$template" | sed "s/{$placeholder}/$value/g")
  printf -v "${out}_path" '%s' "$path"
  printf -v "${out}_directory" '%s' "$(dirname "$path")"
  printf -v "${out}_name" '%s' "$(basename "$path")"
}

# ── Id subsetting ────────────────────────────────────────────────────────────
# Filter TSV $2 in place to its header plus rows whose IID (column 2) appears in
# the id list $1.
subset_to_ids() {
  local ids=$1 tsv=$2
  awk 'FNR==NR{arr[$1];next} (FNR==1 || (FNR>1 && $2 in arr))' "$ids" "$tsv" > "${tsv}.tmp"
  mv "${tsv}.tmp" "$tsv"
}

# ── Covariate value stratification ───────────────────────────────────────────
# Filter the covariate ($3) and phenotype ($4) TSVs in place to the samples whose
# column named $1 holds the value $2 in the covariate file. A no-op when either
# $1 or $2 is empty, i.e. for an unstratified run.
#
# Fails if the column is missing or the covariate stratum is empty
subset_to_covar_value() {
  local col=$1 val=$2 cov=$3 pheno=$4 n
  if [[ -z "$col" || -z "$val" ]]; then
    return 0
  fi

  if ! awk -v col="$col" -v val="$val" '
    FNR==1 {
      for (i = 1; i <= NF; i++) if ($i == col) c = i
      if (!c) {
        print "FATAL: covariate column " col " not found in " FILENAME > "/dev/stderr"
        exit 1
      }
      next
    }
    $c == val { print $2 }' "$cov" > covar_stratum_ids; then
    return 1
  fi

  n=$(awk 'END{print NR}' covar_stratum_ids)
  if [[ "$n" -eq 0 ]]; then
    echo "FATAL: no samples with ${col}=${val} in ${cov}" >&2
    return 1
  fi
  echo "Stratifying on ${col}=${val}: ${n} samples retained"

  subset_to_ids covar_stratum_ids "$cov"
  subset_to_ids covar_stratum_ids "$pheno"
}

# ── Output directory base ────────────────────────────────────────────────────
# Echo the run's base output dir (callers append their own /meta_analysis/... ).
#   $1=type  $2=invnorm  $3=run_prefix
wgs_output_base() {
  if [ "$1" == "Continuous" ]; then
    case "$2" in
      RINT) echo "${3}_rint_wgs" ;;
      MCC)  echo "${3}_mcc_wgs" ;;
      *)    echo "${3}_raw_wgs" ;;
    esac
  else
    echo "${3}_binary_wgs"
  fi
}

# ── Shared FID/IID pairs ─────────────────────────────────────────────────────
# Print the FID/IID pairs present in both the covariate ($1) and phenotype ($2)
# files, on stdout.
build_subset_ids() {
  awk 'FNR==NR {arr[$1"\t"$2]; next} $1"\t"$2 in arr {print $1"\t"$2}' "$1" "$2"
}

# ── Case FID/IID pairs ───────────────────────────────────────────────────────
# Print the FID/IID pairs of the cases for binary phenotype $2 in regenie
# phenotype file $1, on stdout. regenie codes controls 0 and cases 1 (the
# pipeline never passes --1), and missing values as NA.
build_case_ids() {
  local pheno_file=$1 pheno_name=$2
  awk -v pheno="$pheno_name" '
    FNR==1 { for (i = 3; i <= NF; i++) if ($i == pheno) col = i; next }
    col && $col == 1 { print $1"\t"$2 }' "$pheno_file"
}

# ── Control FID/IID pairs ────────────────────────────────────────────────────
# As build_case_ids, for the controls. Samples with a missing (NA) phenotype are
# in neither list, so the two are not complements of the full sample.
build_control_ids() {
  local pheno_file=$1 pheno_name=$2
  awk -v pheno="$pheno_name" '
    FNR==1 { for (i = 3; i <= NF; i++) if ($i == pheno) col = i; next }
    col && $col == 0 { print $1"\t"$2 }' "$pheno_file"
}

# ── Separating-covariate filter for conditioning variants ────────────────────
# Drop conditioning variants that carry no information in a binary-trait model.
# An allele confined to one arm perfectly separates on that covariate while
# contributing nothing to the case/control contrast, so removing it costs no
# adjustment but removes a source of non-convergence and inflated estimates.
# Variants below the threshold are still excluded from *testing* by the caller,
# which keys that off the unfiltered list.
#
# The criterion is the smallest cell of the 2x2 allele-by-outcome table, so it
# catches separation in both directions: an allele seen only in controls and one
# seen only in cases break the fit the same way. Only the case side is common in
# practice, because these are ultra-rare variants and cases are the smaller arm.
#
# Applies only to engines that fit a joint model. COJO projects marginal
# estimates through an LD matrix and has no likelihood to separate, so the caller
# must not apply this filter to a manc-cojo run — it would drop covariates that
# engine handles perfectly well.
#
# Rewrites $6 with the retained IDs and writes the dropped ones to ${6}.dropped.
#   $1=pfile prefix  $2=pheno file  $3=pheno name
#   $4=min minor-allele count required in *each* of cases and controls
#   $5=input ID list  $6=output ID list
#   $7=FID/IID list of the sample regenie will fit
filter_separating_conditional_vars() {
  local pfile=$1 pheno_file=$2 pheno_name=$3 min_ac=$4 in_list=$5 out_list=$6 keep_ids=$7

  if [[ -z $keep_ids ]]; then
    echo "FATAL: filter_separating_conditional_vars needs the analysis-sample id list" >&2
    exit 1
  fi

  : > "${out_list}.dropped"
  if (( min_ac <= 0 )); then
    cp "$in_list" "$out_list"
    return 0
  fi
  if [[ ! -s $in_list ]]; then
    cp "$in_list" "$out_list"
    return 0
  fi

  build_case_ids "$pheno_file" "$pheno_name" > case_ids_all_for_cond_filter
  build_control_ids "$pheno_file" "$pheno_name" > control_ids_all_for_cond_filter
  # Count within the sample regenie actually fits, not the whole phenotype file.
  # The phenotype file is already restricted to the stratum and covariate stratum, but
  # it still carries samples that drop out for want of covariates or step1
  # predictions. A variant whose only case carrier is one of those would
  # otherwise show a count of 1 here and still separate the fit.
  local arm
  for arm in case control; do
    awk 'FNR==NR { keep[$1"\t"$2]; next } ($1"\t"$2) in keep' \
      "$keep_ids" ${arm}_ids_all_for_cond_filter > ${arm}_ids_for_cond_filter
  done
  if [[ ! -s case_ids_for_cond_filter || ! -s control_ids_for_cond_filter ]]; then
    # An arm did not resolve (phenotype column absent, or no 0s/1s in this
    # stratum). Leave the list alone rather than silently emptying it, but say
    # so — a silent no-op here would look identical to a filter that found
    # nothing.
    echo "WARNING: cases or controls missing for '${pheno_name}' in ${pheno_file}; skipping separation filter" >&2
    cp "$in_list" "$out_list"
    return 0
  fi

  # col 5/6 REF_CT/ALT1_CT, matching the --freq counts idiom used elsewhere
  for arm in case control; do
    plink2 --pfile "$pfile" \
      --extract "$in_list" \
      --keep ${arm}_ids_for_cond_filter \
      --freq counts cols=+reffreq,+alt1freq,-maybeprovref \
      --nonfounders \
      --out cond_${arm}_counts
  done

  # A zero anywhere in the allele-by-outcome table separates the fit, so the
  # threshold applies to the smallest of the four cells rather than to the case
  # side alone.
  awk -v min_ac="$min_ac" -v dropped="${out_list}.dropped" '
    FNR == NR { if (FNR > 1) { ref_ctl[$2] = $5; alt_ctl[$2] = $6 } next }
    FNR > 1 {
      if (!($2 in alt_ctl)) {
        print "FATAL: " $2 " missing from control allele counts" > "/dev/stderr"; exit 1
      }
      m = $5
      if ($6 < m) m = $6
      if (ref_ctl[$2] < m) m = ref_ctl[$2]
      if (alt_ctl[$2] < m) m = alt_ctl[$2]
      if (m >= min_ac) print $2; else print $2 > dropped
    }' cond_control_counts.acount cond_case_counts.acount > "$out_list"

  # Every input ID came from the pgen, so anything missing from the .acount is a
  # column-layout surprise rather than an absent variant — fail loudly.
  local n_in n_kept n_dropped
  n_in=$(wc -l < "$in_list")
  n_kept=$(wc -l < "$out_list")
  n_dropped=$(wc -l < "${out_list}.dropped")
  if (( n_kept + n_dropped != n_in )); then
    echo "FATAL: separation filter accounted for $((n_kept + n_dropped))/${n_in} conditioning variants" >&2
    exit 1
  fi
}

# ── Manc-COJO summary-statistic input ────────────────────────────────────────
# Convert REGENIE step-2 output into the GCTA-COJO ".ma" format that manc-cojo
# reads: SNP A1 A2 freq b se p N.
#
# COJO derives every conditional estimate from *marginal* statistics, so this
# must be built from an UNCONDITIONAL regenie run over the same samples and
# variants as the genotypes handed to manc-cojo — COJO assumes the LD reference
# is the GWAS sample itself.
#
# Columns are located by header name, not position: regenie's column set varies
# with the options in use (--af-cc adds N_CASES/N_CONTROLS, for instance), so
# fixed offsets would silently read the wrong field.
#
# p is reconstructed as 10^-LOG10P, which underflows to 0 past ~1e-308 — exactly
# the range the conditioning variants sit in. p takes no part in the conditional
# -effect maths, but manc-cojo validates every numeric field on read and discards
# rows it cannot make sense of, so it is floored at 1e-300 rather than left at 0.
#
# --allow-empty returns 0 rows instead of failing, for mask builds where regenie
# can legitimately test nothing (every mask below --minMAC).
#
#   build_mancojo_ma [--allow-empty] <out.ma> <regenie file> [<regenie file> ...]
build_mancojo_ma() {
  local allow_empty=0
  if [[ "${1:-}" == "--allow-empty" ]]; then allow_empty=1; shift; fi
  local out=$1; shift
  [[ $# -gt 0 ]] || { echo "FATAL: build_mancojo_ma needs at least one regenie file" >&2; exit 1; }

  printf 'SNP A1 A2 freq b se p N\n' > "$out"
  local f
  for f in "$@"; do
    # piping the case directly avoids word-splitting a command out of a variable
    { case "$f" in
        *.gz) gzip -dc "$f" ;;
        *)    cat "$f" ;;
      esac
    } | awk -v OFS=' ' '
      # locate columns from the header of each file
      BEGIN { split("ID ALLELE0 ALLELE1 A1FREQ N BETA SE LOG10P", n_arr, " "); for (i in n_arr) need[n_arr[i]] = 1 }
      /^##/ { next }
      FNR == 1 || $1 == "CHROM" {
        delete col
        for (i = 1; i <= NF; i++) col[$i] = i
        for (k in need) if (!(k in col)) {
          print "FATAL: regenie output missing column " k > "/dev/stderr"; exit 1
        }
        next
      }
      {
        if (col["EXTRA"] && $col["EXTRA"] == "TEST_FAIL") next
        beta = $col["BETA"]; se = $col["SE"]; freq = $col["A1FREQ"]; n = $col["N"]
        log10p = $col["LOG10P"]
        if (se == "NA" || beta == "NA" || freq == "NA" || log10p == "NA" || se + 0 <= 0) next
        # manc-cojo drops freq at the 0/1 boundary; skip rather than emit rows it will reject
        if (freq + 0 <= 0 || freq + 0 >= 1) next
        p = exp(-log10p * log(10))
        if (p <= 1e-300) p = 1e-300
        print $col["ID"], $col["ALLELE1"], $col["ALLELE0"], freq, beta, se, p, n
      }
    ' >> "$out" || { echo "FATAL: failed converting $f to .ma" >&2; exit 1; }
  done

  local n_rows
  n_rows=$(( $(wc -l < "$out") - 1 ))
  if (( n_rows == 0 && ! allow_empty )); then
    echo "FATAL: no usable rows written to $out" >&2
    exit 1
  fi
  echo "$n_rows"
}

# Keep the regenie rows whose A1FREQ is at most <max>, in place, removing the file
# when no data rows remain. Prints the number of rows kept.
#   filter_regenie_max_a1freq <file.regenie.gz> <max>
filter_regenie_max_a1freq() {
  local f=$1 max=$2 n_kept
  gzip -dc "$f" | awk -v max="$max" '
    /^##/ { print; next }
    !fc {
      for (i = 1; i <= NF; i++) if ($i == "A1FREQ") fc = i
      if (!fc) { print "FATAL: no A1FREQ column in regenie output" > "/dev/stderr"; exit 1 }
      print; next
    }
    $fc != "NA" && $fc + 0 <= max + 0' | gzip > "${f}.tmp"
  mv "${f}.tmp" "$f"
  n_kept=$(gzip -dc "$f" | awk '!/^##/ { n++ } END { print n - 1 }')
  if (( n_kept <= 0 )); then
    rm -f "$f"
    n_kept=0
  fi
  echo "$n_kept"
}

# Write .ma rows present in the pvar in manc-cojo's single-cohort .cma.cojo
# layout, with the conditional columns equal to the marginal ones. Stands in for
# manc-cojo when nothing is left to condition on: it errors on an empty
# --cojo-cond list and writes nothing when none of the IDs are in the .ma.
#   write_uncond_cma <ma> <pvar> <chr> <out.cma.cojo>
write_uncond_cma() {
  local ma=$1 pvar=$2 chr=$3 out=$4
  # manc-cojo writes 23 for chrX and upper-cases alleles
  awk -v chr="${chr}" 'BEGIN {
      OFS = "\t"; c = (chr == "X") ? 23 : chr
      print "Chr", "SNP", "bp", "A1", "A2", "freq", "b", "se", "p", "n", "bC", "bC_se", "pC"
    }
    NR == FNR { if ($1 !~ /^#/) pos[$3] = $2; next }
    FNR > 1 && ($1 in pos) {
      print c, $1, pos[$1], toupper($2), toupper($3), $4, $5, $6, $7, $8, $5, $6, $7
    }' "${pvar}" "${ma}" > "${out}"
}

# ── plink2 with nothing left to keep ─────────────────────────────────────────
# Run plink2 with the given arguments, whose --out is <out prefix>. Returns 1 when
# plink2 stopped because no variants were left after its filters, and exits the
# script on any other plink2 error.
#   plink2_allow_no_variants <out prefix> <plink2 args...> || <handle no variants>
plink2_allow_no_variants() {
  local out=$1 rc=0
  shift
  plink2 "$@" || rc=$?
  (( rc == 0 )) && return 0
  grep -q "No variants remaining" "${out}.log" && return 1
  exit "$rc"
}

# ── DNAnexus upload replacing earlier copies ─────────────────────────────────
# Upload <file> into <folder> of <project>, first removing any file of the same
# name there: DNAnexus keeps both, so a retried job would otherwise leave a stale
# copy next to the new one.
#   dx_upload_replace <file> <project> <folder>
dx_upload_replace() {
  local file=$1 project=$2 folder=$3 old
  for old in $(dx find data --path "${project}:/${folder}/" --name "$(basename "$file")" --class file --norecurse --brief); do
    dx rm "${old}"
  done
  dx upload "$file" --path "${project}:/${folder}/"
}

# ── Step1 predictions ────────────────────────────────────────────────────────
# Print the step1 output glob of each comma-separated prefix in $2, rooted at
# $1, one per line, for the caller to fetch.
#   $1=root  $2=prefixes  $3=type  $4=invnorm  $5=stratum
step1_sources() {
  local root=$1 prefixes=$2 type=$3 invnorm=$4 stratum=$5 prefix
  local -a prefix_list
  IFS=',' read -ra prefix_list <<< "$prefixes"
  for prefix in "${prefix_list[@]}"; do
    echo "${root}/${prefix}_${type}_${invnorm}_Step1_${stratum}_*"
  done
}

# Concatenate the fetched pred lists of each prefix in $1 into $6, then fail
# unless every phenotype in $5 has a prediction: regenie skips a phenotype
# missing from the list with only a warning.
#   $1=prefixes  $2=type  $3=invnorm  $4=stratum  $5=phenos  $6=out
merge_step1_pred_lists() {
  local prefixes=$1 type=$2 invnorm=$3 stratum=$4 phenos=$5 out=$6 prefix missing
  local -a prefix_list
  : > "$out"
  IFS=',' read -ra prefix_list <<< "$prefixes"
  for prefix in "${prefix_list[@]}"; do
    cat "${prefix}_${type}_${invnorm}_Step1_${stratum}_pred.list" >> "$out"
  done
  missing=$(tr ',' '\n' <<< "$phenos" | awk 'FNR==NR {have[$1]; next} NF && !($1 in have)' "$out" -)
  if [[ -n "$missing" ]]; then
    echo "FATAL: no step1 predictions for: $(tr '\n' ' ' <<< "$missing")" >&2
    return 1
  fi
}

fetch_step1_sample_file_gcs() {
  local directory=$1 name=$2
  gcloud storage cp "${directory}/${name}.fam" . 2>/dev/null \
    || gcloud storage cp "${directory}/${name}.psam" .
}

fetch_step1_sample_file_dx() {
  local directory=$1 name=$2
  dx download "${directory}/${name}.fam" 2>/dev/null \
    || dx download "${directory}/${name}.psam"
}

strip_fam_header() {
  local prefix=$1
  [[ -f ${prefix}.fam ]] || return 0
  awk 'NR==1 && /FID/ { next } { print }' "${prefix}.fam" > "${prefix}.fam.tmp"
  mv "${prefix}.fam.tmp" "${prefix}.fam"
}

# Print one IID per line, on stdout, from the step-1 sample file: .fam (PLINK1,
# IID in column 2) or .psam (PLINK2, IID column named in the header).
extract_stratum_ids() {
  local prefix=$1
  if [[ -f ${prefix}.fam ]]; then
    awk 'NR==1 && /FID/ { next }   # skip optional header
         { print $2 }' "${prefix}.fam"
  elif [[ -f ${prefix}.psam ]]; then
    # A .psam may carry several header lines; only the last starts with #FID or
    # #IID and it is the one that names the columns. Take the last such line.
    awk '/^#/ {
           if ($1 == "#FID" || $1 == "#IID") {
             col = 0
             for (i = 1; i <= NF; i++) if ($i == "IID" || $i == "#IID") col = i
           }
           header = 1
           next
         }
         !seen {
           seen = 1
           if (!col) {
             if (header) { print "FATAL: no IID column in psam header" > "/dev/stderr"; exit 1 }
             col = 2   # headerless .psam is .fam-shaped
           }
         }
         { print $col }' "${prefix}.psam"
  else
    echo "FATAL: neither ${prefix}.fam nor ${prefix}.psam found" >&2
    exit 1
  fi
}
