//! ACAT of every unit, mask and frequency bin, and conversion of annotations to the
//! parquet schema.
//!
//! `compute_acat_streaming_rs` has two backends:
//! - `"sweep"` (default, fastest): reads `annotations.parquet` in batches, sorted by
//!   unit start, and processes each (unit, chromosome) group in parallel once complete.
//! - `"arrow"`: loads all annotations into memory, so needs no sort order, and at most
//!   64 annotations and 64 masks.
//!
//! `convert_tsv_to_parquet_rs` converts REGENIE annotation and mask files to the schema.

use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use flate2::write::GzEncoder;
use flate2::Compression as GzCompression;

// All backends: csv crate + arrow-parquet + rayon
use arrow_array::{Array, Int64Array, StringArray, RecordBatch};
use arrow_schema::{DataType, Field, Schema};
use flate2::read::GzDecoder;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use rayon::prelude::*;

use crate::acat_core::{acat_log10p_core, weight_standard};

// ── Chromosome order ──────────────────────────────────────────────────

/// Sort rank of a chromosome: 1-22, X (or 23), Y, MT, then others.
#[allow(dead_code)]
fn chrom_rank(c: &str) -> u32 {
    match c {
        "X" | "23" => 23,
        "Y" => 24,
        "MT" => 25,
        s => s.parse::<u32>().unwrap_or(26),
    }
}

// ── Structs ───────────────────────────────────────────────────────────

/// A variant's single-variant result in each frequency bin.
struct VariantResult {
    /// Result per bin, indexed as the thresholds; `None` where the variant is not in the bin.
    per_bin: Vec<Option<BinResult>>,
}

/// A variant's contribution to a bin's ACAT.
struct BinResult {
    /// -log10(P).
    log10p: f64,
    /// Unnormalised weight.
    weight: f64,
}

impl VariantResult {
    /// The variant's result in a bin, if it is in it.
    fn bin(&self, bin_idx: usize) -> Option<&BinResult> {
        self.per_bin.get(bin_idx).and_then(|b| b.as_ref())
    }
}

/// A burden result's contribution to its mask's ACAT.
struct BurdenEntry {
    /// -log10(P).
    log10p: f64,
    /// Unnormalised weight.
    weight: f64,
}

/// An aggregate frequency bin.
struct AafBin {
    /// Upper allele frequency, `label` parsed.
    max_af: f64,
    /// The bin as written for REGENIE `--aaf-bins`, in burden IDs and in the `THRESHOLD` column.
    label: String,
}

impl AafBin {
    /// Parse an aggregate-frequency label (e.g. `"0.001"`) into a bin.
    fn parse(label: &str) -> Result<AafBin, String> {
        label
            .parse::<f64>()
            .map(|max_af| AafBin { max_af, label: label.to_string() })
            .map_err(|_| format!("invalid aggregate_freq '{label}': not a valid number"))
    }
}

/// Column names and separator of the meta-analysis files ACAT reads.
pub struct MetaColumnConfig {
    /// Chromosome column.
    pub chrom_col: String,
    /// ID column.
    pub id_col: String,
    /// Allele frequency column.
    pub a1freq_col: String,
    /// Test column.
    pub test_col: String,
    /// -log10(P) column.
    pub log10p_col: String,
    /// Field separator.
    pub separator: u8,
}

impl Default for MetaColumnConfig {
    /// fedcat's meta-analysis output columns, tab-separated.
    fn default() -> Self {
        Self {
            chrom_col: "CHROM".into(),
            id_col: "ID".into(),
            a1freq_col: "A1FREQ".into(),
            test_col: "TEST".into(),
            log10p_col: "LOG10P".into(),
            separator: b'\t',
        }
    }
}

impl MetaColumnConfig {
    /// Build from a Python dict keyed `chrom`, `id`, `a1freq`, `test`, `log10p` and
    /// `separator`, with defaults for keys left out.
    pub fn from_pydict(dict: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let mut cfg = Self::default();
        if let Some(d) = dict {
            if let Some(v) = d.get_item("chrom")? { cfg.chrom_col = v.extract()?; }
            if let Some(v) = d.get_item("id")? { cfg.id_col = v.extract()?; }
            if let Some(v) = d.get_item("a1freq")? { cfg.a1freq_col = v.extract()?; }
            if let Some(v) = d.get_item("test")? { cfg.test_col = v.extract()?; }
            if let Some(v) = d.get_item("log10p")? { cfg.log10p_col = v.extract()?; }
            if let Some(v) = d.get_item("separator")? {
                let s: String = v.extract()?;
                if let Some(c) = s.bytes().next() {
                    cfg.separator = c;
                }
            }
        }
        Ok(cfg)
    }
}

/// An ACAT output row.
struct AcatResult {
    /// Chromosome.
    chr: String,
    /// First annotated position of the unit's mask; `None` (written NA) without one.
    start: Option<i64>,
    /// Last annotated position of the unit's mask.
    end: Option<i64>,
    /// `UNIT.MASK.BIN`.
    id: String,
    unit: String,
    mask: String,
    /// Frequency bin label.
    threshold: String,
    /// Burden results combined (0 or 1).
    n_agg_combined: usize,
    /// Single-variant results combined.
    n_singlevar_combined: usize,
    /// -log10(P) of the ACAT.
    log10p_acat: f64,
}

// ── CSV/TSV file readers ───────────────────────────────────────────────

/// Read REGENIE mask files (`MASK<tab>annot1,annot2,...`), merging masks of the same name.
///
/// # Returns
/// Mask names in first-seen order, and the annotations of each.
///
/// # Panics
/// If a mask name contains a `.`.
fn load_masks_tsv(paths: Vec<String>) -> (Vec<String>, Vec<Vec<String>>) {
    use std::collections::hash_map::Entry;

    // Annotations per mask, with masks in first-seen order
    let mut mask_order: Vec<String> = Vec::new();
    let mut mask_annots: HashMap<String, Vec<String>> = HashMap::new();

    for path in &paths {
        let reader = open_tsv_reader(path);
        for line in reader.lines() {
            let line = match line {
                Ok(l) => l,
                Err(_) => continue,
            };
            if line.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = line.splitn(2, '\t').collect();
            if parts.len() < 2 {
                continue;
            }
            let name = parts[0].to_string();
            if name.contains('.') {
                panic!("mask name '{name}' contains a dot, which is not allowed");
            }
            let file_annots: Vec<String> = parts[1]
                .split(',')
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
                .collect();

            match mask_annots.entry(name.clone()) {
                Entry::Vacant(e) => {
                    mask_order.push(name);
                    e.insert(file_annots);
                }
                Entry::Occupied(mut e) => {
                    let existing = e.get_mut();
                    for ann in file_annots {
                        if !existing.contains(&ann) {
                            existing.push(ann);
                        }
                    }
                }
            }
        }
    }

    let annots: Vec<Vec<String>> = mask_order.iter().map(|n| mask_annots.remove(n).unwrap()).collect();
    (mask_order, annots)
}

// ═══════════════════════════════════════════════════════════════════════
// Reading inputs, and ACAT per (unit, chromosome)
// ═══════════════════════════════════════════════════════════════════════

/// Read a string column of a record batch.
///
/// # Panics
/// If the column is missing or not a string column.
fn get_string_values(batch: &RecordBatch, col_name: &str) -> Vec<String> {
    let col = batch
        .column_by_name(col_name)
        .unwrap_or_else(|| panic!("missing column {col_name}"));
    if let Some(arr) = col.as_any().downcast_ref::<StringArray>() {
        (0..arr.len()).map(|i| arr.value(i).to_string()).collect()
    } else if let Some(arr) = col
        .as_any()
        .downcast_ref::<arrow_array::LargeStringArray>()
    {
        (0..arr.len()).map(|i| arr.value(i).to_string()).collect()
    } else {
        panic!(
            "column {col_name} is not String or LargeString: {:?}",
            col.data_type()
        );
    }
}

/// Read an integer column of a record batch as i64.
///
/// # Panics
/// If the column is missing or not Int32 or Int64.
fn get_i64_values(batch: &RecordBatch, col_name: &str) -> Vec<i64> {
    let col = batch
        .column_by_name(col_name)
        .unwrap_or_else(|| panic!("missing column {col_name}"));
    if let Some(arr) = col.as_any().downcast_ref::<Int64Array>() {
        (0..arr.len()).map(|i| arr.value(i)).collect()
    } else if let Some(arr) = col.as_any().downcast_ref::<arrow_array::Int32Array>() {
        (0..arr.len()).map(|i| arr.value(i) as i64).collect()
    } else {
        panic!(
            "column {col_name} is not Int32 or Int64: {:?}",
            col.data_type()
        );
    }
}

/// Read single-variant results, keyed by variant code, for the variants in the schema.
///
/// Without an `AAF_BIN` column a variant enters every bin its A1FREQ falls in. With
/// one (a meta-analysis with `aaf_bins`), each row is that bin's result; rows of bins
/// not requested are ignored. Rows with a missing or negative -log10(P) are skipped.
///
/// # Returns
/// Results per variant code, and the set of codes read.
///
/// # Panics
/// If the file cannot be read or lacks a column.
fn load_singles_direct(
    path: &str,
    thresholds: &[AafBin],
    weight_fn: fn(f64) -> f64,
    use_weights: bool,
    col_cfg: &MetaColumnConfig,
    vid_to_code: &HashMap<&str, i64>,
) -> (HashMap<i64, VariantResult>, HashSet<i64>) {
    let reader = open_tsv_reader(path);

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(col_cfg.separator)
        .has_headers(true)
        .comment(Some(b'#'))
        .from_reader(reader);

    let headers = rdr.headers().expect("missing header").clone();
    let find_col = |name: &str| -> usize {
        headers
            .iter()
            .position(|h| h == name)
            .unwrap_or_else(|| {
                panic!(
                    "column '{name}' not found in header: {:?}",
                    headers.iter().collect::<Vec<_>>()
                )
            })
    };
    let id_idx = find_col(&col_cfg.id_col);
    let af_idx = find_col(&col_cfg.a1freq_col);
    let log10p_idx = find_col(&col_cfg.log10p_col);
    let bin_idx_col = headers.iter().position(|h| h == "AAF_BIN");
    let max_threshold = thresholds.iter().map(|b| b.max_af).fold(0.0f64, f64::max);

    let mut code_map: HashMap<i64, VariantResult> = HashMap::new();
    let mut needed: HashSet<i64> = HashSet::new();

    for result in rdr.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };
        let id = match record.get(id_idx) {
            Some(s) => s,
            None => continue,
        };
        let code = match vid_to_code.get(id) {
            Some(&c) => c,
            None => continue,
        };
        let af: f64 = match record.get(af_idx).and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        let log10p: f64 = match record.get(log10p_idx).and_then(|s| s.parse().ok()) {
            Some(v) if v >= 0.0 => v,
            _ => continue,
        };
        let w = if use_weights { weight_fn(af) } else { 1.0 };
        let entry = |vr: &mut VariantResult, i: usize| {
            vr.per_bin[i] = Some(BinResult { log10p, weight: w });
        };
        match bin_idx_col {
            Some(bi) => {
                let bin_af: f64 = match record.get(bi).and_then(|s| s.parse().ok()) {
                    Some(v) => v,
                    None => continue,
                };
                // Match on value, as labels may be written differently ("1e-3", "0.001")
                let Some(i) = thresholds.iter().position(|b| b.max_af == bin_af) else { continue };
                let vr = code_map.entry(code).or_insert_with(|| VariantResult {
                    per_bin: (0..thresholds.len()).map(|_| None).collect(),
                });
                entry(vr, i);
            }
            None => {
                if af > max_threshold {
                    continue;
                }
                let mut vr = VariantResult { per_bin: (0..thresholds.len()).map(|_| None).collect() };
                for (i, b) in thresholds.iter().enumerate() {
                    if af <= b.max_af {
                        entry(&mut vr, i);
                    }
                }
                code_map.insert(code, vr);
            }
        }
        needed.insert(code);
    }
    (code_map, needed)
}

/// Read burden results, keyed by (chromosome, unit, mask, bin label, test).
///
/// IDs are parsed as `UNIT.MASK.BIN` for the given bin labels; other rows are skipped.
///
/// # Panics
/// If the file cannot be read or lacks a column.
fn load_burden_csv(
    path: &str,
    weight_fn: fn(f64) -> f64,
    use_weights: bool,
    thresholds: &[String],
    col_cfg: &MetaColumnConfig,
) -> HashMap<(String, String, String, String, String), BurdenEntry> {
    let reader = open_tsv_reader(path);

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(col_cfg.separator)
        .has_headers(true)
        .comment(Some(b'#'))
        .from_reader(reader);

    let headers = rdr.headers().expect("missing header").clone();
    let find_col = |name: &str| -> usize {
        headers
            .iter()
            .position(|h| h == name)
            .unwrap_or_else(|| {
                panic!(
                    "column '{name}' not found in header: {:?}",
                    headers.iter().collect::<Vec<_>>()
                )
            })
    };
    let id_idx = find_col(&col_cfg.id_col);
    let af_idx = find_col(&col_cfg.a1freq_col);
    let log10p_idx = find_col(&col_cfg.log10p_col);
    let test_idx = find_col(&col_cfg.test_col);
    let chrom_idx = find_col(&col_cfg.chrom_col);

    // UNIT.MASK.BIN, for the bin labels given
    let parse_umt = |id: &str| -> Option<(String, String, String)> {
        for t in thresholds {
            let suffix = format!(".{t}");
            if let Some(prefix) = id.strip_suffix(suffix.as_str()) {
                if let Some(dot_pos) = prefix.rfind('.') {
                    return Some((
                        prefix[..dot_pos].to_string(),
                        prefix[dot_pos + 1..].to_string(),
                        t.clone(),
                    ));
                }
            }
        }
        None
    };

    let mut map = HashMap::new();
    for result in rdr.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };
        let chr = match record.get(chrom_idx) {
            Some(s) => s.to_string(),
            None => continue,
        };
        let burden_id = match record.get(id_idx) {
            Some(s) => s,
            None => continue,
        };
        let (unit, mask, threshold) = match parse_umt(burden_id) {
            Some(t) => t,
            None => continue,
        };
        let test = match record.get(test_idx) {
            Some(s) => s,
            None => continue,
        };
        let af: f64 = match record.get(af_idx).and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        let log10p: f64 = match record.get(log10p_idx).and_then(|s| s.parse().ok()) {
            Some(v) if v >= 0.0 => v,
            _ => continue,
        };
        let w = if use_weights { weight_fn(af) } else { 1.0 };
        map.insert(
            (chr, unit, mask, threshold, test.to_string()),
            BurdenEntry { log10p, weight: w },
        );
    }
    map
}

// ── Parquet schema readers ───────────────────────────────────────────

/// Read every record batch of a parquet file.
///
/// # Panics
/// If the file cannot be read.
fn read_parquet_arrow(path: &str) -> Vec<RecordBatch> {
    let file = File::open(path).unwrap_or_else(|e| panic!("cannot open parquet {path}: {e}"));
    let builder =
        ParquetRecordBatchReaderBuilder::try_new(file).expect("parquet reader build");
    let reader = builder.build().expect("parquet reader");
    reader
        .collect::<Result<Vec<_>, _>>()
        .expect("reading parquet")
}

/// Read `dim_variants.parquet` as parallel code and variant ID lists.
fn load_dim_variants_arrow(schema_dir: &str) -> (Vec<i64>, Vec<String>) {
    let path = format!("{schema_dir}/dim_variants.parquet");
    let batches = read_parquet_arrow(&path);
    let mut codes = Vec::new();
    let mut ids = Vec::new();
    for batch in &batches {
        let c = get_i64_values(batch, "VARIANT_CODE");
        let v = get_string_values(batch, "VARIANT_ID");
        for i in 0..batch.num_rows() {
            codes.push(c[i]);
            ids.push(v[i].clone());
        }
    }
    (codes, ids)
}

/// Read `dim_units.parquet` as parallel code and unit name lists.
fn load_dim_units_arrow(schema_dir: &str) -> (Vec<i64>, Vec<String>) {
    let path = format!("{schema_dir}/dim_units.parquet");
    let batches = read_parquet_arrow(&path);
    let mut codes = Vec::new();
    let mut names = Vec::new();
    for batch in &batches {
        let c = get_i64_values(batch, "UNIT_CODE");
        let n = get_string_values(batch, "UNIT");
        for i in 0..batch.num_rows() {
            codes.push(c[i]);
            names.push(n[i].clone());
        }
    }
    (codes, names)
}

/// Read `dim_annots.parquet` as parallel code and annotation name lists.
fn load_dim_annots_arrow(schema_dir: &str) -> (Vec<i64>, Vec<String>) {
    let path = format!("{schema_dir}/dim_annots.parquet");
    let batches = read_parquet_arrow(&path);
    let mut codes = Vec::new();
    let mut names = Vec::new();
    for batch in &batches {
        let c = get_i64_values(batch, "ANNOT_CODE");
        let n = get_string_values(batch, "ANNOT");
        for i in 0..batch.num_rows() {
            codes.push(c[i]);
            names.push(n[i].clone());
        }
    }
    (codes, names)
}

/// Read `mask_annotations.parquet` as mask names in first-seen order, with each one's annotations.
fn load_mask_annotations_parquet(schema_dir: &str) -> (Vec<String>, Vec<Vec<String>>) {
    let path = format!("{schema_dir}/mask_annotations.parquet");
    let batches = read_parquet_arrow(&path);
    let mut mask_order: Vec<String> = Vec::new();
    let mut mask_annots_map: HashMap<String, Vec<String>> = HashMap::new();
    for batch in &batches {
        let mask_strs = get_string_values(batch, "MASK");
        let annot_strs = get_string_values(batch, "ANNOT");
        for i in 0..batch.num_rows() {
            let mask = mask_strs[i].clone();
            let annot = annot_strs[i].clone();
            if !mask_annots_map.contains_key(&mask) {
                mask_order.push(mask.clone());
                mask_annots_map.insert(mask.clone(), Vec::new());
            }
            mask_annots_map.get_mut(&mask).unwrap().push(annot);
        }
    }
    let annots_vec: Vec<Vec<String>> = mask_order
        .iter()
        .map(|m| mask_annots_map[m].clone())
        .collect();
    (mask_order, annots_vec)
}

/// Read `unit_mask_bounds.parquet`, keyed by (chromosome, unit, mask).
fn load_unit_mask_bounds_arrow(
    schema_dir: &str,
) -> HashMap<(String, String, String), (i64, i64)> {
    let path = format!("{schema_dir}/unit_mask_bounds.parquet");
    let batches = read_parquet_arrow(&path);
    let mut map = HashMap::new();
    for batch in &batches {
        let unit_strs = get_string_values(batch, "UNIT");
        let chr_strs = get_string_values(batch, "CHR");
        let mask_strs = get_string_values(batch, "MASK");
        let min_col = get_i64_values(batch, "MIN_POS");
        let max_col = get_i64_values(batch, "MAX_POS");
        for i in 0..batch.num_rows() {
            let key = (
                chr_strs[i].clone(),
                unit_strs[i].clone(),
                mask_strs[i].clone(),
            );
            map.insert(key, (min_col[i], max_col[i]));
        }
    }
    map
}

/// Read `annotations.parquet` into (unit, chromosome) groups, for the needed variants.
///
/// # Returns
/// `(unit_code, chr, annotations)` per group, where `annotations` maps each variant
/// code to a bitmask of its annotation indices.
fn load_annotation_groups_arrow(
    schema_dir: &str,
    needed_variant_codes: &HashSet<i64>,
    code_to_annot_idx: &HashMap<i64, usize>,
) -> Vec<(i64, String, HashMap<i64, u64>)> {
    let path = format!("{schema_dir}/annotations.parquet");
    let batches = read_parquet_arrow(&path);

    // Chromosome names, interned
    let mut chr_intern: HashMap<String, u16> = HashMap::new();
    let mut chr_names: Vec<String> = Vec::new();

    let mut groups_map: HashMap<(i64, u16), HashMap<i64, u64>> = HashMap::new();
    for batch in &batches {
        let vc = get_i64_values(batch, "VARIANT_CODE");
        let uc = get_i64_values(batch, "UNIT_CODE");
        let ac = get_i64_values(batch, "ANNOT_CODE");
        let chr_strs = get_string_values(batch, "CHR");
        for i in 0..batch.num_rows() {
            let v = vc[i];
            if !needed_variant_codes.contains(&v) {
                continue;
            }
            if let Some(&ai) = code_to_annot_idx.get(&ac[i]) {
                let chr_idx = *chr_intern.entry(chr_strs[i].clone()).or_insert_with(|| {
                    let idx = chr_names.len() as u16;
                    chr_names.push(chr_strs[i].clone());
                    idx
                });
                *groups_map
                    .entry((uc[i], chr_idx))
                    .or_default()
                    .entry(v)
                    .or_insert(0u64) |= 1u64 << ai;
            }
        }
    }

    groups_map
        .into_iter()
        .map(|((uc, ci), va)| (uc, chr_names[ci as usize].clone(), va))
        .collect()
}

/// Compute the ACAT results of one (unit, chromosome) group, for every mask and bin.
///
/// Each mask combines its annotated variants' results in the bin with its burden
/// result, if any. Annotations are sets of indices; used by the sweep backend when
/// there are more than 64 annotations or masks.
///
/// # Arguments
/// * `unit`, `chr` - The group.
/// * `thresholds` - Frequency bins.
/// * `test` - Test of the burden results to combine.
/// * `var_annots` - Annotation indices of each variant in the group.
/// * `variant_results` - Single-variant results by variant code.
/// * `burden_map` - Burden results.
/// * `annot_to_mask_indices` - Masks each annotation index belongs to.
/// * `mask_names` - Mask names, by index.
/// * `unit_mask_bounds` - Positions of each (chromosome, unit, mask).
fn process_unit_chr(
    unit: &str,
    chr: &str,
    thresholds: &[AafBin],
    test: &str,
    var_annots: &HashMap<i64, HashSet<usize>>,
    variant_results: &HashMap<i64, &VariantResult>,
    burden_map: &HashMap<(String, String, String, String, String), BurdenEntry>,
    annot_to_mask_indices: &HashMap<usize, Vec<usize>>,
    mask_names: &[String],
    unit_mask_bounds: &HashMap<(String, String, String), (i64, i64)>,
) -> Vec<AcatResult> {
    let n_masks = mask_names.len();
    let chr_s = chr.to_string();
    let unit_s = unit.to_string();
    let test_s = test.to_string();
    let mut results = Vec::new();

    for (bin_idx, bin) in thresholds.iter().enumerate() {
        let mut mask_log10ps: Vec<Vec<f64>> = vec![Vec::new(); n_masks];
        let mut mask_weights: Vec<Vec<f64>> = vec![Vec::new(); n_masks];

        for (&var_code, annot_set) in var_annots {
            let vr = match variant_results.get(&var_code).and_then(|v| v.bin(bin_idx)) {
                Some(v) => v,
                None => continue,
            };
            if vr.log10p.is_nan() {
                continue;
            }
            let mut hit_masks: HashSet<usize> = HashSet::new();
            for &annot_idx in annot_set {
                if let Some(mask_idxs) = annot_to_mask_indices.get(&annot_idx) {
                    for &mi in mask_idxs {
                        hit_masks.insert(mi);
                    }
                }
            }
            for mi in hit_masks {
                mask_log10ps[mi].push(vr.log10p);
                mask_weights[mi].push(vr.weight);
            }
        }

        let mut burden_masks: HashSet<usize> = HashSet::new();
        for (mi, mn) in mask_names.iter().enumerate() {
            if burden_map.contains_key(&(chr_s.clone(), unit_s.clone(), mn.clone(), bin.label.clone(), test_s.clone())) {
                burden_masks.insert(mi);
            }
        }

        for mi in 0..n_masks {
            let has_vars = !mask_log10ps[mi].is_empty();
            let has_burden = burden_masks.contains(&mi);
            if !has_vars && !has_burden {
                continue;
            }

            let mut lps = std::mem::take(&mut mask_log10ps[mi]);
            let mut wts = std::mem::take(&mut mask_weights[mi]);

            let mut added_burden = false;
            if has_burden {
                let burden_key = (chr_s.clone(), unit_s.clone(), mask_names[mi].clone(), bin.label.clone(), test_s.clone());
                if let Some(be) = burden_map.get(&burden_key) {
                    lps.push(be.log10p);
                    wts.push(be.weight);
                    added_burden = true;
                }
            }

            let wt_sum: f64 = wts.iter().sum();
            if wt_sum <= 0.0 {
                continue;
            }
            for w in &mut wts {
                *w /= wt_sum;
            }

            let log10p_acat = acat_log10p_core(&lps, &wts);
            let acat_id = format!("{}.{}.{}", unit, mask_names[mi], bin.label);

            let (start, end) = unit_mask_bounds
                .get(&(chr_s.clone(), unit_s.clone(), mask_names[mi].clone()))
                .map_or((None, None), |&(s, e)| (Some(s), Some(e)));

            results.push(AcatResult {
                chr: chr_s.clone(),
                start,
                end,
                id: acat_id,
                unit: unit_s.clone(),
                mask: mask_names[mi].clone(),
                threshold: bin.label.clone(),
                n_agg_combined: if added_burden { 1 } else { 0 },
                n_singlevar_combined: lps.len() - if added_burden { 1 } else { 0 },
                log10p_acat,
            });
        }
    }
    results
}

/// As [`process_unit_chr`], with annotations and masks as `u64` bitmasks (at most 64 of each).
///
/// `var_annots` maps each variant to a bitmask of its annotation indices, and
/// `annot_to_mask_bits` gives each annotation index's bitmask of masks.
fn process_unit_chr_bitmask(
    unit: &str,
    chr: &str,
    thresholds: &[AafBin],
    test: &str,
    var_annots: &HashMap<i64, u64>,
    variant_results: &HashMap<i64, &VariantResult>,
    burden_map: &HashMap<(String, String, String, String, String), BurdenEntry>,
    annot_to_mask_bits: &[u64],
    mask_names: &[String],
    unit_mask_bounds: &HashMap<(String, String, String), (i64, i64)>,
) -> Vec<AcatResult> {
    let n_masks = mask_names.len();
    let chr_s = chr.to_string();
    let unit_s = unit.to_string();
    let test_s = test.to_string();
    let mut results = Vec::new();

    for (bin_idx, bin) in thresholds.iter().enumerate() {
        let mut mask_log10ps: Vec<Vec<f64>> = vec![Vec::new(); n_masks];
        let mut mask_weights: Vec<Vec<f64>> = vec![Vec::new(); n_masks];

        for (&var_code, &annot_bits) in var_annots {
            let vr = match variant_results.get(&var_code).and_then(|v| v.bin(bin_idx)) {
                Some(v) => v,
                None => continue,
            };
            if vr.log10p.is_nan() {
                continue;
            }
            // Masks of the variant's annotations
            let mut hit_mask_bits: u64 = 0;
            let mut bits = annot_bits;
            while bits != 0 {
                let ai = bits.trailing_zeros() as usize;
                if ai < annot_to_mask_bits.len() {
                    hit_mask_bits |= annot_to_mask_bits[ai];
                }
                bits &= bits - 1;
            }
            let mut mbits = hit_mask_bits;
            while mbits != 0 {
                let mi = mbits.trailing_zeros() as usize;
                mask_log10ps[mi].push(vr.log10p);
                mask_weights[mi].push(vr.weight);
                mbits &= mbits - 1;
            }
        }

        // Masks with a burden result in this bin
        let mut burden_mask_bits: u64 = 0;
        for (mi, mn) in mask_names.iter().enumerate() {
            if mi < 64
                && burden_map.contains_key(&(
                    chr_s.clone(),
                    unit_s.clone(),
                    mn.clone(),
                    bin.label.clone(),
                    test_s.clone(),
                ))
            {
                burden_mask_bits |= 1u64 << mi;
            }
        }

        for mi in 0..n_masks {
            let has_vars = !mask_log10ps[mi].is_empty();
            let has_burden = burden_mask_bits & (1u64 << mi) != 0;
            if !has_vars && !has_burden {
                continue;
            }

            let mut lps = std::mem::take(&mut mask_log10ps[mi]);
            let mut wts = std::mem::take(&mut mask_weights[mi]);

            let mut added_burden = false;
            if has_burden {
                let burden_key = (
                    chr_s.clone(),
                    unit_s.clone(),
                    mask_names[mi].clone(),
                    bin.label.clone(),
                    test_s.clone(),
                );
                if let Some(be) = burden_map.get(&burden_key) {
                    lps.push(be.log10p);
                    wts.push(be.weight);
                    added_burden = true;
                }
            }

            let wt_sum: f64 = wts.iter().sum();
            if wt_sum <= 0.0 {
                continue;
            }
            for w in &mut wts {
                *w /= wt_sum;
            }

            let log10p_acat = acat_log10p_core(&lps, &wts);
            let acat_id = format!("{}.{}.{}", unit, mask_names[mi], bin.label);

            let (start, end) = unit_mask_bounds
                .get(&(chr_s.clone(), unit_s.clone(), mask_names[mi].clone()))
                .map_or((None, None), |&(s, e)| (Some(s), Some(e)));

            results.push(AcatResult {
                chr: chr_s.clone(),
                start,
                end,
                id: acat_id,
                unit: unit_s.clone(),
                mask: mask_names[mi].clone(),
                threshold: bin.label.clone(),
                n_agg_combined: if added_burden { 1 } else { 0 },
                n_singlevar_combined: lps.len() - if added_burden { 1 } else { 0 },
                log10p_acat,
            });
        }
    }
    results
}

// ═══════════════════════════════════════════════════════════════════════
// Backends
// ═══════════════════════════════════════════════════════════════════════

/// Arrow backend: load all annotation groups into memory and process them in parallel.
///
/// Needs no sort order. Units with burden results but no annotated single-variant
/// results are processed after the groups.
///
/// # Returns
/// `(results, n_singles, n_burden, n_annotation_rows)`.
///
/// # Panics
/// If there are more than 64 annotations or masks, or an input cannot be read.
fn run_arrow_backend(
    singles_tsv: &str,
    burden_tsv: &str,
    schema_dir: &str,
    thresholds: &[AafBin],
    weight_fn: fn(f64) -> f64,
    use_weights: bool,
    col_cfg: &MetaColumnConfig,
) -> (Vec<AcatResult>, usize, usize, usize) {
    let burden_thresholds: Vec<String> = thresholds.iter().map(|b| b.label.clone()).collect();

    // Dimension tables first, to key the single-variant results by code
    let (mask_names, mask_annots) = load_mask_annotations_parquet(schema_dir);
    let (var_codes, var_ids) = load_dim_variants_arrow(schema_dir);
    let (unit_codes, unit_names) = load_dim_units_arrow(schema_dir);
    let (annot_codes, annot_names) = load_dim_annots_arrow(schema_dir);
    // Masks and annotations are bits of a u64 here
    if annot_names.len() > 64 || mask_names.len() > 64 {
        panic!(
            "the arrow ACAT backend handles at most 64 annotations and 64 masks (got {} and {}); \
             use backend='sweep'",
            annot_names.len(),
            mask_names.len()
        );
    }

    let mut vid_to_code: HashMap<&str, i64> = HashMap::with_capacity(var_codes.len());
    for (i, &code) in var_codes.iter().enumerate() {
        vid_to_code.insert(&var_ids[i], code);
    }

    let (code_to_vr_owned, needed_codes) = load_singles_direct(
        singles_tsv, thresholds, weight_fn, use_weights, col_cfg, &vid_to_code,
    );
    let n_singles = code_to_vr_owned.len();
    drop(vid_to_code);

    let burden: HashMap<(String, String, String, String, String), BurdenEntry> =
        if burden_tsv.is_empty() {
            HashMap::new()
        } else {
            load_burden_csv(burden_tsv, weight_fn, use_weights, &burden_thresholds, col_cfg)
        };
    let n_burden = burden.len();

    let mut code_to_unit: HashMap<i64, String> = HashMap::with_capacity(unit_codes.len());
    for (i, &code) in unit_codes.iter().enumerate() {
        code_to_unit.insert(code, unit_names[i].clone());
    }

    let mut code_to_annot_idx: HashMap<i64, usize> = HashMap::with_capacity(annot_codes.len());
    for (i, &code) in annot_codes.iter().enumerate() {
        code_to_annot_idx.insert(code, i);
    }

    // Bitmask of masks per annotation index
    let n_annots_dim = annot_names.len();
    let mut annot_to_mask_bits: Vec<u64> = vec![0u64; n_annots_dim];
    for (mi, mask_ann) in mask_annots.iter().enumerate() {
        for ann in mask_ann {
            if let Some(ai) = annot_names.iter().position(|a| a == ann) {
                annot_to_mask_bits[ai] |= 1u64 << mi;
            }
        }
    }

    let unit_mask_bounds = load_unit_mask_bounds_arrow(schema_dir);

    let groups = load_annotation_groups_arrow(schema_dir, &needed_codes, &code_to_annot_idx);
    let n_annots = groups
        .iter()
        .map(|(_, _, va)| va.values().map(|b| b.count_ones() as usize).sum::<usize>())
        .sum::<usize>();

    // (chr, unit, bin, test) of every group with annotations, to find burden-only units
    let test = "ADD";
    let mut seen_unit_chrs: HashSet<(String, String, String, String)> = HashSet::new();
    for (uc, chr, va) in &groups {
        if va.is_empty() {
            continue;
        }
        if let Some(u) = code_to_unit.get(uc) {
            for bin in thresholds {
                seen_unit_chrs.insert((chr.clone(), u.clone(), bin.label.clone(), test.to_string()));
            }
        }
    }

    let code_to_vr_ref: HashMap<i64, &VariantResult> = code_to_vr_owned
        .iter()
        .map(|(&k, v)| (k, v))
        .collect();

    let mut all_results: Vec<AcatResult> = groups
        .par_iter()
        .filter_map(|(uc, chr, var_annots)| {
            if var_annots.is_empty() {
                return None;
            }
            let unit_name = code_to_unit.get(uc)?;
            Some(process_unit_chr_bitmask(
                unit_name,
                chr,
                thresholds,
                test,
                var_annots,
                &code_to_vr_ref,
                &burden,
                &annot_to_mask_bits,
                &mask_names,
                &unit_mask_bounds,
            ))
        })
        .flatten()
        .collect();

    // Units with burden results but no annotation group, each processed once for every bin
    let mut burden_only: HashSet<(String, String, String)> = HashSet::new();
    for (chr, unit, _mask, threshold, test_key) in burden.keys() {
        let key = (chr.clone(), unit.clone(), threshold.clone(), test_key.clone());
        if !seen_unit_chrs.contains(&key) {
            burden_only.insert((chr.clone(), unit.clone(), test_key.clone()));
        }
    }
    for (chr, unit, test_key) in &burden_only {
        let empty_va: HashMap<i64, u64> = HashMap::new();
        all_results.extend(process_unit_chr_bitmask(
            unit,
            chr,
            thresholds,
            test_key,
            &empty_va,
            &code_to_vr_ref,
            &burden,
            &annot_to_mask_bits,
            &mask_names,
            &unit_mask_bounds,
        ));
    }

    (all_results, n_singles, n_burden, n_annots)
}

/// Process completed bitmask groups in parallel, appending their results and recording them as seen.
fn stream_flush_bitmask(
    completed: &mut Vec<(i64, String, HashMap<i64, u64>)>,
    all_results: &mut Vec<AcatResult>,
    seen: &mut HashSet<(String, String, String, String)>,
    code_to_unit: &HashMap<i64, String>,
    code_to_vr: &HashMap<i64, &VariantResult>,
    burden: &HashMap<(String, String, String, String, String), BurdenEntry>,
    annot_to_mask_bits: &[u64],
    mask_names: &[String],
    unit_mask_bounds: &HashMap<(String, String, String), (i64, i64)>,
    thresholds: &[AafBin],
    test: &str,
) {
    for (uc, chr, va) in completed.iter() {
        if !va.is_empty() {
            if let Some(u) = code_to_unit.get(uc) {
                for bin in thresholds {
                    seen.insert((chr.clone(), u.clone(), bin.label.clone(), test.to_string()));
                }
            }
        }
    }

    let batch_results: Vec<AcatResult> = completed
        .par_iter()
        .filter_map(|(uc, chr, var_annots)| {
            if var_annots.is_empty() {
                return None;
            }
            let unit_name = code_to_unit.get(uc)?;
            Some(process_unit_chr_bitmask(
                unit_name, chr, thresholds, test,
                var_annots, code_to_vr, burden,
                annot_to_mask_bits, mask_names,
                unit_mask_bounds,
            ))
        })
        .flatten()
        .collect();
    all_results.extend(batch_results);
    completed.clear();
}

/// As [`stream_flush_bitmask`], for groups with annotations as index sets.
fn stream_flush_hashset(
    completed: &mut Vec<(i64, String, HashMap<i64, HashSet<usize>>)>,
    all_results: &mut Vec<AcatResult>,
    seen: &mut HashSet<(String, String, String, String)>,
    code_to_unit: &HashMap<i64, String>,
    code_to_vr: &HashMap<i64, &VariantResult>,
    burden: &HashMap<(String, String, String, String, String), BurdenEntry>,
    annot_to_mask_indices: &HashMap<usize, Vec<usize>>,
    mask_names: &[String],
    unit_mask_bounds: &HashMap<(String, String, String), (i64, i64)>,
    thresholds: &[AafBin],
    test: &str,
) {
    for (uc, chr, va) in completed.iter() {
        if !va.is_empty() {
            if let Some(u) = code_to_unit.get(uc) {
                for bin in thresholds {
                    seen.insert((chr.clone(), u.clone(), bin.label.clone(), test.to_string()));
                }
            }
        }
    }

    let batch_results: Vec<AcatResult> = completed
        .par_iter()
        .filter_map(|(uc, chr, var_annots)| {
            if var_annots.is_empty() {
                return None;
            }
            let unit_name = code_to_unit.get(uc)?;
            Some(process_unit_chr(
                unit_name, chr, thresholds, test,
                var_annots, code_to_vr, burden,
                annot_to_mask_indices, mask_names,
                unit_mask_bounds,
            ))
        })
        .flatten()
        .collect();
    all_results.extend(batch_results);
    completed.clear();
}

/// Record a completed (unit, chromosome) group.
///
/// # Panics
/// If it completed before, as then the annotations are not sorted by unit, as the
/// sweep backend needs.
fn check_group_unseen(
    finished_groups: &mut HashSet<(i64, String)>,
    unit_code: i64,
    chr: &str,
    code_to_unit: &HashMap<i64, String>,
) {
    if !finished_groups.insert((unit_code, chr.to_string())) {
        let uname = code_to_unit.get(&unit_code).map(|s| s.as_str()).unwrap_or("?");
        panic!(
            "annotations.parquet is not correctly sorted for sweep backend: \
             unit '{uname}' (code {unit_code}), chr '{chr}' appeared again after being completed. \
             Re-convert with sort_order='unit_start' or use backend='arrow'."
        );
    }
}

/// Sweep backend: stream `annotations.parquet`, processing groups in parallel as they complete.
///
/// The annotations must be sorted by chromosome, unit start, unit and position, as
/// `convert_tsv_to_parquet_rs` writes them, so each (unit, chromosome) group is
/// contiguous. Completed groups are processed in batches of 512. Annotations are
/// bitmasks with at most 64 annotations and masks, else index sets. Units with
/// burden results but no annotated single-variant results are processed last.
///
/// # Returns
/// `(results, n_singles, n_burden, n_annotation_rows)`.
///
/// # Panics
/// If the annotations are not sorted by unit, or an input cannot be read.
fn run_sweep_backend(
    singles_tsv: &str,
    burden_tsv: &str,
    schema_dir: &str,
    thresholds: &[AafBin],
    weight_fn: fn(f64) -> f64,
    use_weights: bool,
    col_cfg: &MetaColumnConfig,
) -> (Vec<AcatResult>, usize, usize, usize) {
    let burden_thresholds: Vec<String> = thresholds.iter().map(|b| b.label.clone()).collect();

    // Dimension tables first, to key the single-variant results by code
    let (mask_names, mask_annots) = load_mask_annotations_parquet(schema_dir);
    let (var_codes, var_ids) = load_dim_variants_arrow(schema_dir);
    let (unit_codes, unit_names) = load_dim_units_arrow(schema_dir);
    let (annot_codes, annot_names) = load_dim_annots_arrow(schema_dir);

    let n_annots_dim = annot_names.len();
    let n_masks = mask_names.len();
    let use_bitmask = n_annots_dim <= 64 && n_masks <= 64;

    let mut vid_to_code: HashMap<&str, i64> = HashMap::with_capacity(var_codes.len());
    for (i, &code) in var_codes.iter().enumerate() {
        vid_to_code.insert(&var_ids[i], code);
    }

    let (code_to_vr_owned, needed_codes) = load_singles_direct(
        singles_tsv, thresholds, weight_fn, use_weights, col_cfg, &vid_to_code,
    );
    let n_singles = code_to_vr_owned.len();
    drop(vid_to_code);

    let burden: HashMap<(String, String, String, String, String), BurdenEntry> =
        if burden_tsv.is_empty() {
            HashMap::new()
        } else {
            load_burden_csv(burden_tsv, weight_fn, use_weights, &burden_thresholds, col_cfg)
        };
    let n_burden = burden.len();

    let code_to_vr_ref: HashMap<i64, &VariantResult> =
        code_to_vr_owned.iter().map(|(&k, v)| (k, v)).collect();

    let mut code_to_unit: HashMap<i64, String> = HashMap::with_capacity(unit_codes.len());
    for (i, &code) in unit_codes.iter().enumerate() {
        code_to_unit.insert(code, unit_names[i].clone());
    }

    let mut code_to_annot_idx: HashMap<i64, usize> = HashMap::with_capacity(annot_codes.len());
    for (i, &code) in annot_codes.iter().enumerate() {
        code_to_annot_idx.insert(code, i);
    }

    // Bitmask of masks per annotation index, with at most 64 of each
    let annot_to_mask_bits: Vec<u64> = if use_bitmask {
        let mut bits = vec![0u64; n_annots_dim];
        for (mi, mask_ann) in mask_annots.iter().enumerate() {
            for ann in mask_ann {
                if let Some(ai) = annot_names.iter().position(|a| a == ann) {
                    bits[ai] |= 1u64 << mi;
                }
            }
        }
        bits
    } else {
        Vec::new()
    };

    // Otherwise mask indices per annotation index
    let annot_to_mask_indices: HashMap<usize, Vec<usize>> = if !use_bitmask {
        let mut map: HashMap<usize, Vec<usize>> = HashMap::new();
        for (mi, mask_ann) in mask_annots.iter().enumerate() {
            for ann in mask_ann {
                if let Some(ai) = annot_names.iter().position(|a| a == ann) {
                    map.entry(ai).or_default().push(mi);
                }
            }
        }
        map
    } else {
        HashMap::new()
    };

    let unit_mask_bounds = load_unit_mask_bounds_arrow(schema_dir);

    // ── Stream annotations.parquet ──
    let annot_path = format!("{schema_dir}/annotations.parquet");
    let annot_file =
        File::open(&annot_path).unwrap_or_else(|e| panic!("cannot open parquet {annot_path}: {e}"));
    let annot_reader = ParquetRecordBatchReaderBuilder::try_new(annot_file)
        .expect("parquet reader build")
        .build()
        .expect("parquet reader");

    let test = "ADD";
    let mut all_results: Vec<AcatResult> = Vec::new();
    let mut seen_unit_chrs: HashSet<(String, String, String, String)> = HashSet::new();
    let mut n_annots: usize = 0;

    // The group being read, and completed groups waiting to be processed
    let mut pending_bm: Option<(i64, String, HashMap<i64, u64>)> = None;
    let mut completed_bm: Vec<(i64, String, HashMap<i64, u64>)> = Vec::new();
    let mut pending_hs: Option<(i64, String, HashMap<i64, HashSet<usize>>)> = None;
    let mut completed_hs: Vec<(i64, String, HashMap<i64, HashSet<usize>>)> = Vec::new();

    const FLUSH_THRESHOLD: usize = 512;

    // Completed groups, to catch a unit split by unsorted annotations
    let mut finished_groups: HashSet<(i64, String)> = HashSet::new();

    for batch_result in annot_reader {
        let batch = batch_result.expect("reading parquet batch");
        let vc = get_i64_values(&batch, "VARIANT_CODE");
        let uc = get_i64_values(&batch, "UNIT_CODE");
        let ac = get_i64_values(&batch, "ANNOT_CODE");
        let chr_strs = get_string_values(&batch, "CHR");

        for i in 0..batch.num_rows() {
            let v = vc[i];
            if !needed_codes.contains(&v) {
                continue;
            }
            let ai = match code_to_annot_idx.get(&ac[i]) {
                Some(&ai) => ai,
                None => continue,
            };
            n_annots += 1;

            let unit_code = uc[i];
            let chr = &chr_strs[i];

            if use_bitmask {
                let is_new = match &pending_bm {
                    Some((pu, pc, _)) => *pu != unit_code || pc != chr,
                    None => true,
                };
                if is_new {
                    if let Some(prev) = pending_bm.take() {
                        check_group_unseen(&mut finished_groups, prev.0, &prev.1, &code_to_unit);
                        completed_bm.push(prev);
                    }
                    pending_bm = Some((unit_code, chr.clone(), HashMap::new()));
                }
                if let Some((_, _, ref mut va)) = pending_bm {
                    *va.entry(v).or_insert(0u64) |= 1u64 << ai;
                }
            } else {
                let is_new = match &pending_hs {
                    Some((pu, pc, _)) => *pu != unit_code || pc != chr,
                    None => true,
                };
                if is_new {
                    if let Some(prev) = pending_hs.take() {
                        check_group_unseen(&mut finished_groups, prev.0, &prev.1, &code_to_unit);
                        completed_hs.push(prev);
                    }
                    pending_hs = Some((unit_code, chr.clone(), HashMap::new()));
                }
                if let Some((_, _, ref mut va)) = pending_hs {
                    va.entry(v).or_default().insert(ai);
                }
            }
        }

        if use_bitmask && completed_bm.len() >= FLUSH_THRESHOLD {
            stream_flush_bitmask(
                &mut completed_bm,
                &mut all_results,
                &mut seen_unit_chrs,
                &code_to_unit,
                &code_to_vr_ref,
                &burden,
                &annot_to_mask_bits,
                &mask_names,
                &unit_mask_bounds,
                thresholds,
                test,
            );
        } else if !use_bitmask && completed_hs.len() >= FLUSH_THRESHOLD {
            stream_flush_hashset(
                &mut completed_hs,
                &mut all_results,
                &mut seen_unit_chrs,
                &code_to_unit,
                &code_to_vr_ref,
                &burden,
                &annot_to_mask_indices,
                &mask_names,
                &unit_mask_bounds,
                thresholds,
                test,
            );
        }
    }

    if let Some(prev) = pending_bm.take() {
        check_group_unseen(&mut finished_groups, prev.0, &prev.1, &code_to_unit);
        completed_bm.push(prev);
    }
    if let Some(prev) = pending_hs.take() {
        check_group_unseen(&mut finished_groups, prev.0, &prev.1, &code_to_unit);
        completed_hs.push(prev);
    }

    if use_bitmask {
        stream_flush_bitmask(
            &mut completed_bm,
            &mut all_results,
            &mut seen_unit_chrs,
            &code_to_unit,
            &code_to_vr_ref,
            &burden,
            &annot_to_mask_bits,
            &mask_names,
            &unit_mask_bounds,
            thresholds,
            test,
        );
    } else {
        stream_flush_hashset(
            &mut completed_hs,
            &mut all_results,
            &mut seen_unit_chrs,
            &code_to_unit,
            &code_to_vr_ref,
            &burden,
            &annot_to_mask_indices,
            &mask_names,
            &unit_mask_bounds,
            thresholds,
            test,
        );
    }

    // Units with burden results but no annotation group, each processed once for every bin
    let mut burden_only: HashSet<(String, String, String)> = HashSet::new();
    for (chr, unit, _mask, threshold, test_key) in burden.keys() {
        let key = (chr.clone(), unit.clone(), threshold.clone(), test_key.clone());
        if !seen_unit_chrs.contains(&key) {
            burden_only.insert((chr.clone(), unit.clone(), test_key.clone()));
        }
    }
    if use_bitmask {
        for (chr, unit, test_key) in &burden_only {
            let empty_va: HashMap<i64, u64> = HashMap::new();
            all_results.extend(process_unit_chr_bitmask(
                unit, chr, thresholds, test_key, &empty_va, &code_to_vr_ref, &burden,
                &annot_to_mask_bits, &mask_names, &unit_mask_bounds,
            ));
        }
    } else {
        for (chr, unit, test_key) in &burden_only {
            let empty_va: HashMap<i64, HashSet<usize>> = HashMap::new();
            all_results.extend(process_unit_chr(
                unit, chr, thresholds, test_key, &empty_va, &code_to_vr_ref, &burden,
                &annot_to_mask_indices, &mask_names, &unit_mask_bounds,
            ));
        }
    }

    (all_results, n_singles, n_burden, n_annots)
}

// ═══════════════════════════════════════════════════════════════════════
// Python function
// ═══════════════════════════════════════════════════════════════════════

/// Compute the ACAT of every unit, mask and frequency bin, and write them as a TSV.
///
/// Each ACAT combines a mask's single-variant results in the bin with its burden
/// result, if any. Burden results are those of test `ADD`, with IDs `UNIT.MASK.BIN`.
///
/// # Arguments
/// * `singles_tsv` - Single-variant meta-analysis results.
/// * `burden_tsv` - Burden (MAC-10) meta-analysis results, or `""` for none.
/// * `schema_dir` - Annotation schema directory.
/// * `aggregate_freqs` - Frequency bin labels, e.g. `["0.001"]`.
/// * `weighting_scheme` - `"standard"` ([`weight_standard`]) or `"equal"`.
/// * `output_path` - Output TSV, or `""` to write nothing.
/// * `column_config` - Input column names; see [`MetaColumnConfig::from_pydict`].
/// * `compression` - `"zst"`, `"gz"`, or anything else for none.
/// * `compression_level` - Compression level; `None` for 3 (zst) or 6 (gz).
/// * `backend` - `"sweep"`, or `"arrow"` (any other value means sweep).
///
/// # Returns
/// A dict of counts: `_n_singles`, `_n_burden`, `_n_annots` and `_n_results`.
///
/// # Errors
/// `ValueError` for no or invalid frequency bins. Unreadable inputs panic, which
/// Python receives as a `PanicException`.
#[pyfunction]
#[pyo3(signature = (
    singles_tsv,
    burden_tsv,
    schema_dir,
    aggregate_freqs,
    weighting_scheme = "standard",
    output_path = "",
    column_config = None,
    compression = "zst",
    compression_level = None,
    backend = "sweep",
))]
pub fn compute_acat_streaming_rs(
    py: Python<'_>,
    singles_tsv: &str,
    burden_tsv: &str,
    schema_dir: &str,
    aggregate_freqs: Vec<String>,
    weighting_scheme: &str,
    output_path: &str,
    column_config: Option<&Bound<'_, PyDict>>,
    compression: &str,
    compression_level: Option<i32>,
    backend: &str,
) -> PyResult<PyObject> {
    let col_cfg = MetaColumnConfig::from_pydict(column_config)?;
    let weight_fn: fn(f64) -> f64 = match weighting_scheme {
        "standard" => weight_standard,
        _ => |_| 1.0,
    };
    let use_weights = weighting_scheme != "equal";

    if aggregate_freqs.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "aggregate_freqs must contain at least one threshold",
        ));
    }
    let thresholds: Vec<AafBin> = aggregate_freqs
        .iter()
        .map(|s| AafBin::parse(s).map_err(pyo3::exceptions::PyValueError::new_err))
        .collect::<PyResult<Vec<_>>>()?;

    // Release the GIL for the computation and output: callers run this off the event loop
    let backend = backend.to_string();
    let output_path = output_path.to_string();
    let compression = compression.to_string();
    let (n, n_singles, n_burden, n_annots) = py.allow_threads(move || {
    let (all_results, n_singles, n_burden, n_annots) = match backend.as_str() {
        "arrow" => run_arrow_backend(
            singles_tsv,
            burden_tsv,
            schema_dir,
            &thresholds,
            weight_fn,
            use_weights,
            &col_cfg,
        ),
        _ => run_sweep_backend(
            singles_tsv,
            burden_tsv,
            schema_dir,
            &thresholds,
            weight_fn,
            use_weights,
            &col_cfg,
        ),
    };

    let n = all_results.len();

    if !output_path.is_empty() {
        let file = File::create(&output_path)
            .unwrap_or_else(|e| panic!("cannot create output file {output_path}: {e}"));

        macro_rules! write_rows {
            ($writer:expr) => {{
                let w = &mut $writer;
                writeln!(w, "CHR\tSTART\tEND\tID\tUNIT\tMASK\tTHRESHOLD\tN_AGG_COMBINED\tN_SINGLEVAR_COMBINED\tLOG10P_ACAT")
                    .unwrap_or_else(|e| panic!("write header: {e}"));
                for r in &all_results {
                    let p_str = if r.log10p_acat.is_nan() {
                        "NA".to_string()
                    } else {
                        format!("{}", r.log10p_acat)
                    };
                    let pos_str = |p: Option<i64>| p.map_or_else(|| "NA".to_string(), |v| v.to_string());
                    writeln!(
                        w,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        r.chr, pos_str(r.start), pos_str(r.end), r.id, r.unit, r.mask, r.threshold,
                        r.n_agg_combined, r.n_singlevar_combined, p_str
                    )
                    .unwrap_or_else(|e| panic!("write row: {e}"));
                }
            }};
        }

        match compression.as_str() {
            "gz" => {
                let level = compression_level.map(|l| l as u32).unwrap_or(6);
                let enc = GzEncoder::new(file, GzCompression::new(level));
                let mut writer = BufWriter::new(enc);
                write_rows!(writer);
                writer
                    .into_inner()
                    .unwrap_or_else(|e| panic!("flush gz: {e}"))
                    .finish()
                    .unwrap_or_else(|e| panic!("finish gz: {e}"));
            }
            "zst" => {
                let level = compression_level.unwrap_or(3);
                let enc = zstd::Encoder::new(file, level)
                    .unwrap_or_else(|e| panic!("zstd init: {e}"))
                    .auto_finish();
                let mut writer = BufWriter::new(enc);
                write_rows!(writer);
            }
            _ => {
                let mut writer = BufWriter::new(file);
                write_rows!(writer);
            }
        }
    }
    (n, n_singles, n_burden, n_annots)
    });

    let dict = PyDict::new(py);
    dict.set_item("_n_singles", n_singles)?;
    dict.set_item("_n_burden", n_burden)?;
    dict.set_item("_n_annots", n_annots)?;
    dict.set_item("_n_results", n)?;

    Ok(dict.into())
}

// ═══════════════════════════════════════════════════════════════════════
// TSV → Parquet conversion
// ═══════════════════════════════════════════════════════════════════════

/// Open a text file, plain, gzip (`.gz`) or zstd (`.zst`, `.zstd`), for reading lines.
///
/// # Panics
/// If the file cannot be opened.
fn open_tsv_reader(path: &str) -> Box<dyn BufRead> {
    let file = File::open(path).unwrap_or_else(|e| panic!("cannot open TSV {path}: {e}"));
    if path.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else if path.ends_with(".zst") || path.ends_with(".zstd") {
        let dec = zstd::Decoder::new(file)
            .unwrap_or_else(|e| panic!("cannot decode zstd TSV {path}: {e}"));
        Box::new(BufReader::new(dec))
    } else {
        Box::new(BufReader::new(file))
    }
}

/// Sort rank of a chromosome: 1-22, X, Y, MT, then others.
fn chr_rank(chr: &str) -> u32 {
    match chr {
        "X" | "23" => 23,
        "Y" | "24" => 24,
        "MT" | "M" | "25" => 25,
        _ => chr.parse::<u32>().unwrap_or(26),
    }
}

/// Write a dimension table of codes 0.. and names to parquet.
fn write_dim_parquet(path: &str, code_col: &str, name_col: &str, names: &[String]) {
    let codes: Vec<i64> = (0..names.len() as i64).collect();
    let code_arr = Int64Array::from(codes);
    let name_arr = StringArray::from(names.to_vec());
    let schema = Schema::new(vec![
        Field::new(code_col, DataType::Int64, false),
        Field::new(name_col, DataType::Utf8, false),
    ]);
    let batch = RecordBatch::try_new(
        std::sync::Arc::new(schema),
        vec![std::sync::Arc::new(code_arr), std::sync::Arc::new(name_arr)],
    )
    .expect("build dim batch");
    let file = File::create(path).unwrap_or_else(|e| panic!("cannot create {path}: {e}"));
    let mut writer = ArrowWriter::try_new(file, batch.schema(), None).expect("arrow writer");
    writer.write(&batch).expect("write dim batch");
    writer.close().expect("close dim writer");
}

/// Whether an ID is a `chr:pos:ref:alt` variant ID, with non-empty parts and a numeric position.
#[inline]
fn is_valid_vid(s: &str) -> bool {
    let b = s.as_bytes();
    let len = b.len();
    if len < 7 {
        return false; // shortest: "X:1:A:T"
    }
    let mut colons = [0usize; 3];
    let mut n_colons = 0;
    for (i, &ch) in b.iter().enumerate() {
        if ch == b':' {
            if n_colons >= 3 {
                return false;
            }
            colons[n_colons] = i;
            n_colons += 1;
        }
    }
    if n_colons != 3 {
        return false;
    }
    if colons[0] == 0 || colons[1] == colons[0] + 1 || colons[2] == colons[1] + 1 || colons[2] == len - 1 {
        return false;
    }
    for &ch in &b[colons[0] + 1..colons[1]] {
        if !ch.is_ascii_digit() {
            return false;
        }
    }
    true
}

/// Reorder five parallel arrays in place by a permutation, so `a[i] = old_a[perm[i]]`.
///
/// `perm` is left as the identity.
fn apply_permutation_5(
    a1: &mut [i64],
    a2: &mut [i64],
    a3: &mut [i64],
    a4: &mut [u16],
    a5: &mut [i64],
    perm: &mut [usize],
) {
    for idx in 0..perm.len() {
        if perm[idx] == idx {
            continue;
        }
        let mut current = idx;
        let t1 = a1[idx];
        let t2 = a2[idx];
        let t3 = a3[idx];
        let t4 = a4[idx];
        let t5 = a5[idx];
        loop {
            let from = perm[current];
            perm[current] = current;
            if from == idx {
                a1[current] = t1;
                a2[current] = t2;
                a3[current] = t3;
                a4[current] = t4;
                a5[current] = t5;
                break;
            }
            a1[current] = a1[from];
            a2[current] = a2[from];
            a3[current] = a3[from];
            a4[current] = a4[from];
            a5[current] = a5[from];
            current = from;
        }
    }
}

/// Convert REGENIE annotation and mask files to the annotation parquet schema.
///
/// Annotation files are headerless `variant<tab>unit<tab>annotation` (plain, gzip or
/// zstd); rows whose ID is not `chr:pos:ref:alt` are skipped. Codes are assigned in
/// first-seen order, and `annotations.parquet` is sorted by chromosome, unit start,
/// unit and position, as the sweep ACAT backend needs.
///
/// # Arguments
/// * `annotation_tsvs` - REGENIE annotation files.
/// * `masks_tsvs` - REGENIE mask files.
/// * `output_dir` - Schema directory, created if needed.
/// * `valid_variant_ids` - Variant IDs to keep; `None` for all.
/// * `full_bounds` - With `valid_variant_ids`, take mask positions from all annotated
///   variants rather than only those kept.
///
/// # Returns
/// A dict of counts `n_rows`, `n_variants`, `n_units`, `n_annots`, and `full_bounds`,
/// whether all variants' positions were used.
///
/// # Panics
/// If an input cannot be read or an output written.
#[pyfunction]
#[pyo3(signature = (annotation_tsvs, masks_tsvs, output_dir, valid_variant_ids=None, full_bounds=true))]
pub fn convert_tsv_to_parquet_rs(
    py: Python<'_>,
    annotation_tsvs: Vec<String>,
    masks_tsvs: Vec<String>,
    output_dir: &str,
    valid_variant_ids: Option<Vec<String>>,
    full_bounds: bool,
) -> PyResult<PyObject> {
    // Without the GIL, as callers run this in a thread
    let output_dir = output_dir.to_string();
    let (n_rows, n_variants, n_units, n_annots, scan_bounds) = py.allow_threads(move || {
        convert_tsv_to_parquet(annotation_tsvs, masks_tsvs, &output_dir, valid_variant_ids, full_bounds)
    });

    let dict = PyDict::new(py);
    dict.set_item("n_rows", n_rows)?;
    dict.set_item("n_variants", n_variants)?;
    dict.set_item("n_units", n_units)?;
    dict.set_item("n_annots", n_annots)?;
    dict.set_item("full_bounds", scan_bounds)?;

    Ok(dict.into())
}

/// Body of [`convert_tsv_to_parquet_rs`].
///
/// # Returns
/// `(n_rows, n_variants, n_units, n_annots, full_bounds_used)`.
fn convert_tsv_to_parquet(
    annotation_tsvs: Vec<String>,
    masks_tsvs: Vec<String>,
    output_dir: &str,
    valid_variant_ids: Option<Vec<String>>,
    full_bounds: bool,
) -> (usize, usize, usize, usize, bool) {
    std::fs::create_dir_all(output_dir)
        .unwrap_or_else(|e| panic!("cannot create output dir {output_dir}: {e}"));

    let valid_set: Option<HashSet<String>> = valid_variant_ids.map(|ids| ids.into_iter().collect());

    // Codes of units, annotations and variants
    let mut unit_map: HashMap<String, i64> = HashMap::new();
    let mut annot_map: HashMap<String, i64> = HashMap::new();
    let mut variant_map: HashMap<String, i64> = HashMap::new();

    // Chromosome codes, by raw and prefix-stripped name
    let mut chr_map: HashMap<String, u16> = HashMap::new();
    let mut chr_intern: Vec<String> = Vec::new();

    // Rows kept
    let mut unit_codes: Vec<i64> = Vec::new();
    let mut annot_codes: Vec<i64> = Vec::new();
    let mut variant_codes: Vec<i64> = Vec::new();
    let mut chr_codes: Vec<u16> = Vec::new();
    let mut pos_vals: Vec<i64> = Vec::new();

    // With full bounds, position bounds over every row, kept or not, keyed by codes
    // separate from the schema's so dropped units and annotations stay out of it
    let scan_bounds = full_bounds && valid_set.is_some();
    let mut scan_unit_map: HashMap<String, u32> = HashMap::new();
    let mut scan_annot_map: HashMap<String, u32> = HashMap::new();
    let mut scan_minmax: HashMap<(u16, u32, u32), (i64, i64)> = HashMap::new();

    for tsv_path in &annotation_tsvs {
        let reader = open_tsv_reader(tsv_path);
        for line in reader.lines() {
            let line = match line {
                Ok(l) => l,
                Err(_) => continue,
            };
            if line.starts_with('#') {
                continue;
            }
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            // Split on the first two tabs
            let bytes = trimmed.as_bytes();
            let tab1 = match bytes.iter().position(|&b| b == b'\t') {
                Some(p) => p,
                None => continue,
            };
            let tab2 = match bytes[tab1 + 1..].iter().position(|&b| b == b'\t') {
                Some(p) => tab1 + 1 + p,
                None => continue,
            };
            let vid = &trimmed[..tab1];
            let unit_raw = &trimmed[tab1 + 1..tab2];
            let annot = &trimmed[tab2 + 1..];

            if !is_valid_vid(vid) {
                continue;
            }
            let keep = valid_set.as_ref().map_or(true, |vset| vset.contains(vid));
            if !keep && !scan_bounds {
                continue;
            }

            let colon1 = vid.find(':').unwrap();
            let colon2 = vid[colon1 + 1..].find(':').unwrap() + colon1 + 1;
            let chr_part = &vid[..colon1];
            let pos_str = &vid[colon1 + 1..colon2];
            let pos: i64 = pos_str.parse().unwrap_or_else(|_| {
                panic!("invalid position in variant ID: {vid}")
            });

            // Chromosome without a chr prefix, upper case
            let chr_code = if let Some(&code) = chr_map.get(chr_part) {
                code
            } else {
                let stripped = chr_part
                    .strip_prefix("chr")
                    .or_else(|| chr_part.strip_prefix("CHR"))
                    .or_else(|| chr_part.strip_prefix("Chr"))
                    .unwrap_or(chr_part);
                let chr_norm = stripped.to_uppercase();
                let code = chr_intern.len() as u16;
                chr_intern.push(chr_norm.clone());
                chr_map.insert(chr_part.to_string(), code);
                if chr_part != stripped {
                    chr_map.insert(stripped.to_string(), code);
                }
                code
            };

            if scan_bounds {
                let su = if let Some(&c) = scan_unit_map.get(unit_raw) {
                    c
                } else {
                    let next = scan_unit_map.len() as u32;
                    scan_unit_map.insert(unit_raw.to_string(), next);
                    next
                };
                let sa = if let Some(&c) = scan_annot_map.get(annot) {
                    c
                } else {
                    let next = scan_annot_map.len() as u32;
                    scan_annot_map.insert(annot.to_string(), next);
                    next
                };
                scan_minmax
                    .entry((chr_code, su, sa))
                    .and_modify(|(mn, mx)| {
                        *mn = (*mn).min(pos);
                        *mx = (*mx).max(pos);
                    })
                    .or_insert((pos, pos));
            }
            if !keep {
                continue;
            }

            let uc = if let Some(&c) = unit_map.get(unit_raw) {
                c
            } else {
                let next = unit_map.len() as i64;
                unit_map.insert(unit_raw.to_string(), next);
                next
            };
            let ac = if let Some(&c) = annot_map.get(annot) {
                c
            } else {
                let next = annot_map.len() as i64;
                annot_map.insert(annot.to_string(), next);
                next
            };
            let vc = if let Some(&c) = variant_map.get(vid) {
                c
            } else {
                let next = variant_map.len() as i64;
                variant_map.insert(vid.to_string(), next);
                next
            };

            unit_codes.push(uc);
            annot_codes.push(ac);
            variant_codes.push(vc);
            chr_codes.push(chr_code);
            pos_vals.push(pos);
        }
    }

    let n_rows = unit_codes.len();

    // Names by code
    let mut unit_names = vec![String::new(); unit_map.len()];
    for (name, &code) in &unit_map {
        unit_names[code as usize] = name.clone();
    }
    let mut annot_names = vec![String::new(); annot_map.len()];
    for (name, &code) in &annot_map {
        annot_names[code as usize] = name.clone();
    }
    let mut variant_names = vec![String::new(); variant_map.len()];
    for (name, &code) in &variant_map {
        variant_names[code as usize] = name.clone();
    }

    // ── Sort by chromosome, unit start, unit, position ──
    if n_rows > 0 {
        let chr_ranks: Vec<u32> = chr_codes.iter().map(|&c| chr_rank(&chr_intern[c as usize])).collect();

        // Start of each unit on each chromosome
        let mut unit_chr_min: HashMap<(i64, u32), i64> = HashMap::new();
        for i in 0..n_rows {
            let key = (unit_codes[i], chr_ranks[i]);
            let entry = unit_chr_min.entry(key).or_insert(pos_vals[i]);
            if pos_vals[i] < *entry {
                *entry = pos_vals[i];
            }
        }

        let mut indices: Vec<usize> = (0..n_rows).collect();
        indices.sort_unstable_by(|&a, &b| {
            chr_ranks[a]
                .cmp(&chr_ranks[b])
                .then_with(|| {
                    let min_a = unit_chr_min[&(unit_codes[a], chr_ranks[a])];
                    let min_b = unit_chr_min[&(unit_codes[b], chr_ranks[b])];
                    min_a.cmp(&min_b)
                })
                .then_with(|| unit_codes[a].cmp(&unit_codes[b]))
                .then_with(|| pos_vals[a].cmp(&pos_vals[b]))
        });

        apply_permutation_5(
            &mut unit_codes,
            &mut annot_codes,
            &mut variant_codes,
            &mut chr_codes,
            &mut pos_vals,
            &mut indices,
        );
    }

    write_dim_parquet(
        &format!("{output_dir}/dim_variants.parquet"),
        "VARIANT_CODE",
        "VARIANT_ID",
        &variant_names,
    );
    write_dim_parquet(
        &format!("{output_dir}/dim_units.parquet"),
        "UNIT_CODE",
        "UNIT",
        &unit_names,
    );
    write_dim_parquet(
        &format!("{output_dir}/dim_annots.parquet"),
        "ANNOT_CODE",
        "ANNOT",
        &annot_names,
    );

    // ── annotations.parquet ──
    {
        let chr_strs: Vec<&str> = chr_codes.iter().map(|&c| chr_intern[c as usize].as_str()).collect();
        let schema = Schema::new(vec![
            Field::new("VARIANT_CODE", DataType::Int64, false),
            Field::new("UNIT_CODE", DataType::Int64, false),
            Field::new("ANNOT_CODE", DataType::Int64, false),
            Field::new("CHR", DataType::Utf8, false),
            Field::new("POS", DataType::Int64, false),
        ]);
        let batch = RecordBatch::try_new(
            std::sync::Arc::new(schema),
            vec![
                std::sync::Arc::new(Int64Array::from(variant_codes)),
                std::sync::Arc::new(Int64Array::from(unit_codes.clone())),
                std::sync::Arc::new(Int64Array::from(annot_codes.clone())),
                std::sync::Arc::new(StringArray::from(chr_strs)),
                std::sync::Arc::new(Int64Array::from(pos_vals.clone())),
            ],
        )
        .expect("build annotations batch");
        let path = format!("{output_dir}/annotations.parquet");
        let file = File::create(&path).unwrap_or_else(|e| panic!("cannot create {path}: {e}"));
        let mut writer = ArrowWriter::try_new(file, batch.schema(), None).expect("arrow writer");
        writer.write(&batch).expect("write annotations");
        writer.close().expect("close annotations writer");
    }

    // ── unit_mask_bounds.parquet ──
    {
        let (mask_names, mask_annots) = load_masks_tsv(masks_tsvs);

        // Mask indices per annotation code
        let mut annot_to_masks: HashMap<i64, Vec<usize>> = HashMap::new();
        for (mi, mask_ann) in mask_annots.iter().enumerate() {
            for ann in mask_ann {
                if let Some(&ac) = annot_map.get(ann) {
                    annot_to_masks.entry(ac).or_default().push(mi);
                }
            }
        }

        let mut bounds_rows: Vec<(String, String, String, i64, i64)> = Vec::new();
        if scan_bounds {
            let mut scan_annot_to_masks: HashMap<u32, Vec<usize>> = HashMap::new();
            for (mi, mask_ann) in mask_annots.iter().enumerate() {
                for ann in mask_ann {
                    if let Some(&sa) = scan_annot_map.get(ann.as_str()) {
                        scan_annot_to_masks.entry(sa).or_default().push(mi);
                    }
                }
            }
            let mut scan_unit_names = vec![""; scan_unit_map.len()];
            for (name, &code) in &scan_unit_map {
                scan_unit_names[code as usize] = name.as_str();
            }
            let mut mask_bounds: HashMap<(u16, u32, usize), (i64, i64)> = HashMap::new();
            for (&(chr_code, su, sa), &(mn, mx)) in &scan_minmax {
                if let Some(masks) = scan_annot_to_masks.get(&sa) {
                    for &mi in masks {
                        mask_bounds
                            .entry((chr_code, su, mi))
                            .and_modify(|(bmn, bmx)| {
                                *bmn = (*bmn).min(mn);
                                *bmx = (*bmx).max(mx);
                            })
                            .or_insert((mn, mx));
                    }
                }
            }
            for ((chr_code, su, mi), (mn, mx)) in mask_bounds {
                bounds_rows.push((
                    scan_unit_names[su as usize].to_string(),
                    chr_intern[chr_code as usize].clone(),
                    mask_names[mi].clone(),
                    mn,
                    mx,
                ));
            }
        }

        // Otherwise bounds of the kept rows, by (unit, chromosome) group
        let mut i = 0;
        while !scan_bounds && i < n_rows {
            let uc = unit_codes[i];
            let chr = &chr_intern[chr_codes[i] as usize];
            let mut j = i + 1;
            while j < n_rows && unit_codes[j] == uc && chr_intern[chr_codes[j] as usize] == *chr {
                j += 1;
            }
            let mut mask_min: HashMap<usize, i64> = HashMap::new();
            let mut mask_max: HashMap<usize, i64> = HashMap::new();
            for k in i..j {
                let ac = annot_codes[k];
                let pos = pos_vals[k];
                if let Some(masks) = annot_to_masks.get(&ac) {
                    for &mi in masks {
                        mask_min.entry(mi).and_modify(|v| *v = (*v).min(pos)).or_insert(pos);
                        mask_max.entry(mi).and_modify(|v| *v = (*v).max(pos)).or_insert(pos);
                    }
                }
            }
            let unit_name = &unit_names[uc as usize];
            for (&mi, &mn) in &mask_min {
                bounds_rows.push((
                    unit_name.clone(),
                    chr.clone(),
                    mask_names[mi].clone(),
                    mn,
                    mask_max[&mi],
                ));
            }
            i = j;
        }

        let units_col: Vec<&str> = bounds_rows.iter().map(|r| r.0.as_str()).collect();
        let chrs_col: Vec<&str> = bounds_rows.iter().map(|r| r.1.as_str()).collect();
        let masks_col: Vec<&str> = bounds_rows.iter().map(|r| r.2.as_str()).collect();
        let min_col: Vec<i64> = bounds_rows.iter().map(|r| r.3).collect();
        let max_col: Vec<i64> = bounds_rows.iter().map(|r| r.4).collect();

        let schema = Schema::new(vec![
            Field::new("UNIT", DataType::Utf8, false),
            Field::new("CHR", DataType::Utf8, false),
            Field::new("MASK", DataType::Utf8, false),
            Field::new("MIN_POS", DataType::Int64, false),
            Field::new("MAX_POS", DataType::Int64, false),
        ]);
        let batch = RecordBatch::try_new(
            std::sync::Arc::new(schema),
            vec![
                std::sync::Arc::new(StringArray::from(units_col)),
                std::sync::Arc::new(StringArray::from(chrs_col)),
                std::sync::Arc::new(StringArray::from(masks_col)),
                std::sync::Arc::new(Int64Array::from(min_col)),
                std::sync::Arc::new(Int64Array::from(max_col)),
            ],
        )
        .expect("build bounds batch");
        let path = format!("{output_dir}/unit_mask_bounds.parquet");
        let file = File::create(&path).unwrap_or_else(|e| panic!("cannot create {path}: {e}"));
        let mut writer = ArrowWriter::try_new(file, batch.schema(), None).expect("arrow writer");
        writer.write(&batch).expect("write bounds");
        writer.close().expect("close bounds writer");

        // ── mask_annotations.parquet ──
        // A mask with none of its annotations kept gets a placeholder row (ANNOT_CODE
        // -1, which matches no annotation), so ACAT still finds it for burden results
        let mut ma_masks: Vec<String> = Vec::new();
        let mut ma_codes: Vec<i64> = Vec::new();
        let mut ma_names: Vec<String> = Vec::new();
        for (mi, mask_ann) in mask_annots.iter().enumerate() {
            let mut n_found = 0usize;
            for ann in mask_ann {
                if let Some(&ac) = annot_map.get(ann) {
                    ma_masks.push(mask_names[mi].clone());
                    ma_codes.push(ac);
                    ma_names.push(ann.clone());
                    n_found += 1;
                }
            }
            if n_found == 0 {
                ma_masks.push(mask_names[mi].clone());
                ma_codes.push(-1);
                ma_names.push("".to_string());
            }
        }
        {
            let schema = Schema::new(vec![
                Field::new("MASK", DataType::Utf8, false),
                Field::new("ANNOT_CODE", DataType::Int64, false),
                Field::new("ANNOT", DataType::Utf8, false),
            ]);
            let ma_masks_ref: Vec<&str> = ma_masks.iter().map(|s| s.as_str()).collect();
            let ma_names_ref: Vec<&str> = ma_names.iter().map(|s| s.as_str()).collect();
            let batch = RecordBatch::try_new(
                std::sync::Arc::new(schema),
                vec![
                    std::sync::Arc::new(StringArray::from(ma_masks_ref)),
                    std::sync::Arc::new(Int64Array::from(ma_codes)),
                    std::sync::Arc::new(StringArray::from(ma_names_ref)),
                ],
            )
            .expect("build mask_annotations batch");
            let path = format!("{output_dir}/mask_annotations.parquet");
            let file = File::create(&path).unwrap_or_else(|e| panic!("cannot create {path}: {e}"));
            let mut writer = ArrowWriter::try_new(file, batch.schema(), None).expect("arrow writer");
            writer.write(&batch).expect("write mask_annotations");
            writer.close().expect("close mask_annotations writer");
        }
    }

    (n_rows, variant_map.len(), unit_map.len(), annot_map.len(), scan_bounds)
}
