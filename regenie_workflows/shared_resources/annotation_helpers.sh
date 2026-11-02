#!/bin/bash
# Shared annotation-manifest helpers for the aggregate-family regenie scripts
# (step2 aggregates and aggregate conditioning, on both platforms). Source this
# file after downloading the manifest's triplet files; the helpers are not meant
# to be run directly. Platform-specific downloads stay in the callers.

# ── Manifest parsing ─────────────────────────────────────────────────────────
# Parse the annotation_manifest string ("name,mask,annot,setlist;name2,...") into
# the TSV file manifest_resolved, resolving {chr}. Echoes the underscore-joined
# list of triplet names (annots_running) on stdout.
#   annots_running=$(parse_annotation_manifest "$annotation_manifest" "$chr")
parse_annotation_manifest() {
  local manifest=$1 chr=$2
  echo "$manifest" | tr ';' '\n' | tr ',' '\t' | sed "s/{chr}/${chr}/g" > manifest_resolved
  cut -f 1 manifest_resolved | paste -sd "_" -
}

# ── Local-name manifest ──────────────────────────────────────────────────────
# From manifest_resolved (full paths), write manifest_local with the local
# basenames (any .zst suffix stripped) under ./input_files/<name>/, matching
# where the callers download the triplet files.
build_manifest_local() {
  while IFS=$'\t' read -r name mask_path annot_path setlist_path; do
    printf '%s\t%s\t%s\t%s\n' \
      "${name}" \
      "./input_files/${name}/$(basename "$mask_path" .zst)" \
      "./input_files/${name}/$(basename "$annot_path" .zst)" \
      "./input_files/${name}/$(basename "$setlist_path" .zst)"
  done < manifest_resolved > manifest_local
}

# ── Parallel .zst decompression ──────────────────────────────────────────────
# Decompress each passed .zst file (args) to its name minus .zst, in parallel,
# then remove the originals. Returns 1 if any decompression failed so the caller
# can emit its platform-specific error and exit. Assumes plink2 is on PATH.
decompress_zst_parallel() {
  local pids=() fail=0 file out pid
  for file in "$@"; do
    out="${file%.zst}"
    plink2 --zst-decompress "$file" "$out" &
    pids+=("$!")
  done
  for pid in "${pids[@]}"; do
    wait "$pid" || fail=1
  done
  (( fail )) && return 1
  (( $# )) && rm -f -- "$@"
  return 0
}
