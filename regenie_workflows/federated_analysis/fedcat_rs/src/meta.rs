//! Z-score to -log10(P) conversion, and the streaming IVW meta-analysis.
//!
//! - `z_to_log10p_rs`: Z-scores to -log10 of their two-sided p-values, by a Mills
//!   ratio expansion for extreme Z.
//! - `meta_analyse_regenie_streaming_rs`: IVW fixed-effects meta-analysis of REGENIE
//!   results, accumulating one cohort at a time with Polars (about 200-400 MB peak
//!   for 18 cohorts of 10M variants).
//!
//! Also the input column configuration shared with `meta_merge.rs`.

use numpy::{PyArray1, PyReadonlyArray1};
use polars::prelude::*;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use statrs::distribution::{ChiSquared, ContinuousCDF, Normal};

use flate2::write::GzEncoder;
use flate2::Compression as GzCompression;
use std::fs::File;
use std::io::Write;
use std::sync::Arc;


const LN10: f64 = std::f64::consts::LN_10;
/// ln(2π)
const LN_2PI: f64 = 1.837_877_066_409_345_5;
const LN_2: f64 = std::f64::consts::LN_2;

/// Convert a Z-score to -log10 of its two-sided p-value, precise at any |Z|.
pub fn z_to_log10p_scalar(z: f64) -> f64 {
    let az = z.abs();
    if az < 1e-15 {
        return 0.0;
    }

    if az <= 37.0 {
        // sf via erfc keeps full relative precision, which 1 - cdf loses from |Z| ~ 6
        let normal = Normal::new(0.0, 1.0).unwrap();
        let sf = normal.sf(az);
        if sf <= 0.0 {
            return _asymptotic_log10p(az);
        }
        let result = -(sf.ln() + LN_2) / LN10;
        return result.max(0.0);
    }

    // sf underflows beyond 37
    _asymptotic_log10p(az)
}

/// -log10(two-sided P) of a large |z| by the Mills ratio expansion:
/// ln(SF) ≈ -z²/2 - ln(z) - ½ln(2π) + ln(1 - 1/z² + 3/z⁴).
#[inline]
fn _asymptotic_log10p(az: f64) -> f64 {
    let z2 = az * az;
    let inv_z2 = 1.0 / z2;
    let correction = (-inv_z2 + 3.0 * inv_z2 * inv_z2).ln_1p();
    let log_sf = -0.5 * z2 - az.ln() - 0.5 * LN_2PI + correction;
    let result = -(log_sf + LN_2) / LN10;
    result.max(0.0)
}

/// Convert an array of Z-scores to -log10 of their two-sided p-values.
///
/// # Returns
/// A new array of the same length.
#[pyfunction]
pub fn z_to_log10p_rs<'py>(
    py: Python<'py>,
    z_array: PyReadonlyArray1<'py, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let z_slice = z_array.as_slice()?;
    let result: Vec<f64> = z_slice.iter().map(|&z| z_to_log10p_scalar(z)).collect();
    Ok(PyArray1::from_vec(py, result))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_z_zero() {
        assert!((z_to_log10p_scalar(0.0) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_z_196() {
        let result = z_to_log10p_scalar(1.96);
        assert!((result - 1.301).abs() < 0.01);
    }

    #[test]
    fn test_negative_symmetric() {
        let pos = z_to_log10p_scalar(3.0);
        let neg = z_to_log10p_scalar(-3.0);
        assert!((pos - neg).abs() < 1e-10);
    }

    #[test]
    fn test_extreme_z() {
        let result = z_to_log10p_scalar(40.0);
        assert!(result > 300.0);
        assert!(result.is_finite());
    }

    #[test]
    fn test_matches_scipy_logsf() {
        // -(norm.logsf(z) + ln2) / ln10
        for (z, expected) in [
            (6.0, 8.704834331812725),
            (7.9, 14.55454614571694),
            (8.0, 14.905112555353172),
            (20.0, 88.2590653474116),
            (36.9, 297.3352594298938),
        ] {
            let result = z_to_log10p_scalar(z);
            assert!((result - expected).abs() < 1e-9, "z={z}: {result} vs {expected}");
        }
    }

    #[test]
    fn test_very_extreme_z() {
        let result = z_to_log10p_scalar(100.0);
        assert!(result > 2000.0);
        assert!(result.is_finite());
    }
}

// ── Configurable column schema ────────────────────────────────────────────────

/// Column names and separator of meta-analysis input files, defaulting to REGENIE's.
pub struct ColumnConfig {
    /// Chromosome column.
    pub chrom: String,
    /// Position column.
    pub genpos: String,
    /// Variant or aggregate ID column.
    pub id: String,
    /// Reference allele column.
    pub allele0: String,
    /// Tested allele column.
    pub allele1: String,
    /// Tested allele frequency column.
    pub a1freq: String,
    /// Sample size column.
    pub n: String,
    /// Test column.
    pub test: String,
    /// Effect size column.
    pub beta: String,
    /// Standard error column.
    pub se: String,
    /// Column flagging failed tests (`TEST_FAIL`); `None` if there is none.
    pub extra: Option<String>,
    /// -log10(P) column; read from the config but not used.
    pub log10p: Option<String>,
    /// Field separator; `None` to detect tab or space.
    pub separator: Option<char>,
}

impl Default for ColumnConfig {
    /// REGENIE's column names, with the separator detected.
    fn default() -> Self {
        Self {
            chrom: "CHROM".into(),
            genpos: "GENPOS".into(),
            id: "ID".into(),
            allele0: "ALLELE0".into(),
            allele1: "ALLELE1".into(),
            a1freq: "A1FREQ".into(),
            n: "N".into(),
            test: "TEST".into(),
            beta: "BETA".into(),
            se: "SE".into(),
            extra: Some("EXTRA".into()),
            log10p: Some("LOG10P".into()),
            separator: None,
        }
    }
}

impl ColumnConfig {
    /// Build from a Python dict keyed by field name, with defaults for keys left out.
    ///
    /// An empty `extra` or `log10p` sets it to `None`.
    pub fn from_pydict(dict: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let mut cfg = Self::default();
        if let Some(d) = dict {
            if let Some(v) = d.get_item("chrom")? { cfg.chrom = v.extract()?; }
            if let Some(v) = d.get_item("genpos")? { cfg.genpos = v.extract()?; }
            if let Some(v) = d.get_item("id")? { cfg.id = v.extract()?; }
            if let Some(v) = d.get_item("allele0")? { cfg.allele0 = v.extract()?; }
            if let Some(v) = d.get_item("allele1")? { cfg.allele1 = v.extract()?; }
            if let Some(v) = d.get_item("a1freq")? { cfg.a1freq = v.extract()?; }
            if let Some(v) = d.get_item("n")? { cfg.n = v.extract()?; }
            if let Some(v) = d.get_item("test")? { cfg.test = v.extract()?; }
            if let Some(v) = d.get_item("beta")? { cfg.beta = v.extract()?; }
            if let Some(v) = d.get_item("se")? { cfg.se = v.extract()?; }
            if let Some(v) = d.get_item("extra")? {
                let s: String = v.extract()?;
                cfg.extra = if s.is_empty() { None } else { Some(s) };
            }
            if let Some(v) = d.get_item("log10p")? {
                let s: String = v.extract()?;
                cfg.log10p = if s.is_empty() { None } else { Some(s) };
            }
            if let Some(v) = d.get_item("separator")? {
                let s: String = v.extract()?;
                cfg.separator = s.chars().next();
            }
        }
        Ok(cfg)
    }
}

/// Return the configured separator, or detect it from a (possibly compressed) file's header.
///
/// Tab is chosen if it splits out the ID or BETA column, else space.
fn detect_sep_from_scan(path: &str, cfg: &ColumnConfig) -> u8 {
    if let Some(sep) = cfg.separator {
        return sep as u8;
    }
    for &sep in &[b'\t', b' '] {
        if let Ok(mut lf) = LazyCsvReader::new(path.into())
            .with_separator(sep)
            .with_has_header(true)
            .with_comment_prefix(Some("##".into()))
            .with_infer_schema_length(Some(1))
            .with_n_rows(Some(1))
            .with_glob(false)
            .finish()
        {
            if let Ok(schema) = lf.collect_schema() {
                let id: PlSmallStr = cfg.id.as_str().into();
                let beta: PlSmallStr = cfg.beta.as_str().into();
                if schema.contains(&id) || schema.contains(&beta) {
                    return sep;
                }
            }
        }
    }
    b' '
}

// ── Known schema & Polars-based pipeline ─────────────────────────────────────

/// Return the dtypes of the input columns, by their configured names, so Polars infers no types.
fn regenie_dtype_overrides(cfg: &ColumnConfig) -> Arc<Schema> {
    let mut s = Schema::with_capacity(13);
    s.insert(cfg.chrom.as_str().into(), DataType::String);
    s.insert(cfg.genpos.as_str().into(), DataType::Int64);
    s.insert(cfg.id.as_str().into(), DataType::String);
    s.insert(cfg.allele0.as_str().into(), DataType::String);
    s.insert(cfg.allele1.as_str().into(), DataType::String);
    s.insert(cfg.a1freq.as_str().into(), DataType::Float64);
    s.insert(cfg.n.as_str().into(), DataType::Float64);
    s.insert(cfg.test.as_str().into(), DataType::String);
    s.insert(cfg.beta.as_str().into(), DataType::Float64);
    s.insert(cfg.se.as_str().into(), DataType::Float64);
    if let Some(ref extra) = cfg.extra {
        s.insert(extra.as_str().into(), DataType::String);
    }
    Arc::new(s)
}

/// Values read as null.
const REGENIE_NULL_VALUES: [&str; 5] = ["NA", "nan", "-nan", "NaN", "-NaN"];

/// [`REGENIE_NULL_VALUES`] for every column.
fn regenie_null_values() -> NullValues {
    NullValues::AllColumns(REGENIE_NULL_VALUES.iter().map(|s| (*s).into()).collect())
}

/// Scan one REGENIE file, compressed or not, with the given dtypes and no type inference.
fn scan_csv_polars(path: &str, sep: u8, dtypes: &Arc<Schema>) -> PolarsResult<LazyFrame> {
    LazyCsvReader::new(path.into())
        .with_separator(sep)
        .with_has_header(true)
        .with_comment_prefix(Some("##".into()))
        .with_infer_schema_length(Some(0))
        .with_dtype_overwrite(Some(dtypes.clone()))
        .with_missing_is_null(true)
        .with_null_values(Some(regenie_null_values()))
        .with_glob(false)
        .finish()
}

/// Scan one REGENIE file into standard columns, keeping its valid rows.
///
/// Columns are renamed to REGENIE's names. Rows are dropped where SE <= 0, a numeric
/// value is missing or not finite, or EXTRA is `TEST_FAIL`; chromosome 23 becomes X.
fn scan_regenie_df(
    path: &str,
    sep: u8,
    cfg: &ColumnConfig,
    dtypes: &Arc<Schema>,
) -> PolarsResult<LazyFrame> {
    let mut lf = scan_csv_polars(path, sep, dtypes)?;

    let schema = lf.collect_schema()?;
    let extra_src_name = cfg.extra.as_deref().unwrap_or("EXTRA");
    let has_extra = schema.contains(extra_src_name);

    let mut renames_old: Vec<PlSmallStr> = Vec::new();
    let mut renames_new: Vec<PlSmallStr> = Vec::new();
    for (src, dst) in [
        (&cfg.id, "ID"),
        (&cfg.allele0, "ALLELE0"),
        (&cfg.allele1, "ALLELE1"),
        (&cfg.a1freq, "A1FREQ"),
        (&cfg.n, "N"),
        (&cfg.beta, "BETA"),
        (&cfg.se, "SE"),
    ] {
        if src != dst {
            renames_old.push(src.as_str().into());
            renames_new.push(dst.into());
        }
    }
    if cfg.chrom != "CHROM" {
        renames_old.push(cfg.chrom.as_str().into());
        renames_new.push("CHROM".into());
    }
    if cfg.genpos != "GENPOS" {
        renames_old.push(cfg.genpos.as_str().into());
        renames_new.push("GENPOS".into());
    }
    if cfg.test != "TEST" {
        renames_old.push(cfg.test.as_str().into());
        renames_new.push("TEST".into());
    }
    if let Some(ref extra) = cfg.extra {
        if extra != "EXTRA" && has_extra {
            renames_old.push(extra.as_str().into());
            renames_new.push("EXTRA".into());
        }
    }

    let mut lf = lf;

    if !renames_old.is_empty() {
        lf = lf.rename(renames_old, renames_new, true);
    }

    lf = lf.filter(
        col("SE")
            .gt(lit(0.0f64))
            .and(col("SE").is_finite())
            .and(col("A1FREQ").is_finite())
            .and(col("BETA").is_finite())
            .and(col("N").is_finite()),
    );

    if has_extra {
        lf = lf.filter(
            col("EXTRA")
                .neq(lit("TEST_FAIL"))
                .or(col("EXTRA").is_null()),
        );
    }

    lf = lf.with_columns([when(col("CHROM").eq(lit("23")))
        .then(lit("X"))
        .otherwise(col("CHROM"))
        .alias("CHROM")]);

    let select_exprs = vec![
        col("CHROM"),
        col("GENPOS"),
        col("ID"),
        col("TEST"),
        col("ALLELE0"),
        col("ALLELE1"),
        col("A1FREQ"),
        col("N"),
        col("BETA"),
        col("SE"),
    ];
    lf = lf.select(select_exprs);

    Ok(lf)
}

// ── Statistics as Polars expressions ─────────────────────────────────────────
//
// Each wraps a scalar Rust function in one Expr::map node, rather than a deep
// expression tree with many intermediate columns.

/// Output field of a map: the input's name, as Float64.
fn float64_output(_schema: &Schema, field: &Field) -> PolarsResult<Field> {
    Ok(Field::new(field.name().clone(), DataType::Float64))
}

/// `LOG10P` from the `Z` column, by [`z_to_log10p_scalar`].
fn log10p_map_expr() -> Expr {
    col("Z")
        .map(
            |s| {
                let z = s.as_materialized_series().f64()?;
                let out: Float64Chunked = z.apply_values(z_to_log10p_scalar);
                Ok(out.into_column())
            },
            float64_output,
        )
        .alias("LOG10P")
}

/// `Q_P`, the p-value of Cochran's Q on N_STUDIES - 1 degrees of freedom; null for one study.
fn q_p_map_expr() -> Expr {
    as_struct(vec![col("Q"), col("N_STUDIES")])
        .map(
            |s| {
                let ca = s.as_materialized_series().struct_()?;
                let q_f = ca.field_by_name("Q")?.f64()?.clone();
                let ns_s = ca.field_by_name("N_STUDIES")?;
                let ns_cast = ns_s.cast(&DataType::Float64)?;
                let ns_f = ns_cast.f64()?;
                let out: Float64Chunked = q_f
                    .into_iter()
                    .zip(ns_f.into_iter())
                    .map(|(q, ns)| match (q, ns) {
                        (Some(q), Some(ns)) => {
                            let df = ns - 1.0;
                            if df >= 1.0 {
                                Some(ChiSquared::new(df).unwrap().sf(q))
                            } else {
                                None
                            }
                        }
                        _ => None,
                    })
                    .collect();
                Ok(out.into_column())
            },
            float64_output,
        )
        .alias("Q_P")
}

/// `I2`, the I² heterogeneity index (Q - df) / Q clipped to [0, 1]; 0 for one study.
fn i2_map_expr() -> Expr {
    as_struct(vec![col("Q"), col("N_STUDIES")])
        .map(
            |s| {
                let ca = s.as_materialized_series().struct_()?;
                let q_f = ca.field_by_name("Q")?.f64()?.clone();
                let ns_s = ca.field_by_name("N_STUDIES")?;
                let ns_cast = ns_s.cast(&DataType::Float64)?;
                let ns_f = ns_cast.f64()?;
                let out: Float64Chunked = q_f
                    .into_iter()
                    .zip(ns_f.into_iter())
                    .map(|(q, ns)| match (q, ns) {
                        // One study has no heterogeneity; its Q is only rounding error
                        (Some(_), Some(ns)) if ns <= 1.0 => Some(0.0),
                        (Some(q), Some(ns)) => {
                            let df = ns - 1.0;
                            Some(if q > 0.0 {
                                ((q - df) / q).clamp(0.0, 1.0)
                            } else {
                                0.0
                            })
                        }
                        _ => None,
                    })
                    .collect();
                Ok(out.into_column())
            },
            float64_output,
        )
        .alias("I2")
}

/// Largest difference in BETA or SE allowed between a cohort's rows of one variant.
pub const DUPLICATE_TOLERANCE: f64 = 1e-5;

/// Check a cohort's rows of each variant agree, from its per-variant BETA and SE ranges.
///
/// # Errors
/// `ValueError` naming the first variant whose rows differ by more than
/// [`DUPLICATE_TOLERANCE`] in BETA or SE.
fn check_duplicates_agree(cohort_df: &DataFrame, label: &str) -> PyResult<()> {
    let polars_err = |e: PolarsError| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{e}"));
    let conflicts = cohort_df
        .clone()
        .lazy()
        .filter(
            (col("_beta_max") - col("_beta_min"))
                .gt(lit(DUPLICATE_TOLERANCE))
                .or((col("_se_max") - col("_se_min")).gt(lit(DUPLICATE_TOLERANCE))),
        )
        .limit(1)
        .collect()
        .map_err(polars_err)?;
    if conflicts.height() == 0 {
        return Ok(());
    }
    let s = |name: &str| -> PyResult<String> {
        let v = conflicts.column(name).map_err(polars_err)?.get(0).map_err(polars_err)?;
        Ok(v.str_value().to_string())
    };
    Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
        "Conflicting results for variant '{}' (test={}) in cohort '{}': BETA {} to {}, SE {} to {} across its rows",
        s("ID")?, s("TEST")?, label, s("_beta_min")?, s("_beta_max")?, s("_se_min")?, s("_se_max")?,
    )))
}

// ── Allele-frequency bins ────────────────────────────────────────────────────

/// Parse allele-frequency bin labels into `(label, max_af)` pairs.
///
/// With bins, a cohort contributes to bin `t` only where its own A1FREQ is `<= t`, as
/// ACAT within one cohort would, and the output has a row per variant and bin.
///
/// # Errors
/// `ValueError` if a label is not a number in (0, 1].
pub fn parse_aaf_bins(aaf_bins: Option<Vec<String>>) -> PyResult<Vec<(String, f64)>> {
    let mut bins = Vec::new();
    for label in aaf_bins.unwrap_or_default() {
        let max_af: f64 = label.parse().map_err(|_| {
            PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("invalid aaf bin '{label}'"))
        })?;
        if !(max_af > 0.0 && max_af <= 1.0) {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "aaf bin '{label}' must be in (0, 1]"
            )));
        }
        bins.push((label, max_af));
    }
    Ok(bins)
}

// ── Per-cohort Polars accumulation — fast I/O, low memory ────────────────────

/// Meta-analyse REGENIE results by IVW fixed effects, accumulating one cohort at a time.
///
/// Each cohort's files are read and reduced to one row per variant: a variant's rows
/// in several files, or repeated in one, must agree in BETA and SE within
/// [`DUPLICATE_TOLERANCE`], and the first is used. Its IVW terms are added to a
/// running total. BETA, SE, Z, LOG10P, Q, Q_P and I2 are derived from the totals,
/// and rows written sorted by chromosome (1-22, X, Y, MT, others), position, ID and
/// test, after a `## Cohorts:` comment line. A1FREQ is the N-weighted mean, and
/// DIRECTIONS has a `+`, `-` or `?` (absent) per cohort.
///
/// # Arguments
/// * `cohort_files` - `(label, files)` per cohort.
/// * `output_path` - Output TSV.
/// * `res_type` - Result type, for the log only.
/// * `column_config` - Input column names, keyed by [`ColumnConfig`] field.
/// * `compression` - `"zst"`, `"gz"`, or anything else for none.
/// * `compression_level` - Compression level; `None` for 3 (zst) or 6 (gz).
/// * `aaf_bins` - Allele-frequency bins; see [`parse_aaf_bins`].
///
/// # Returns
/// `(n_rows_written, n_cohorts)`.
#[pyfunction]
#[pyo3(signature = (cohort_files, output_path, res_type, column_config=None, compression="zst", compression_level=None, aaf_bins=None))]
pub fn meta_analyse_regenie_streaming_rs(
    py: Python<'_>,
    cohort_files: Vec<(String, Vec<String>)>,
    output_path: &str,
    res_type: &str,
    column_config: Option<&Bound<'_, PyDict>>,
    compression: &str,
    compression_level: Option<i32>,
    aaf_bins: Option<Vec<String>>,
) -> PyResult<(usize, usize)> {
    let cfg = ColumnConfig::from_pydict(column_config)?;
    let bins = parse_aaf_bins(aaf_bins)?;
    let n_cohorts = cohort_files.len();

    let first_file = cohort_files.iter()
        .flat_map(|(_, files)| files.first())
        .next()
        .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyValueError, _>("No input files"))?;

    let t_total = std::time::Instant::now();
    let dtypes = regenie_dtype_overrides(&cfg);
    let sep = detect_sep_from_scan(first_file, &cfg);

    // Variant key, also by bin when binning
    let group_cols = || -> Vec<Expr> {
        let mut g = vec![col("CHROM"), col("GENPOS"), col("ID"), col("TEST")];
        if !bins.is_empty() {
            g.push(col("AAF_BIN"));
        }
        g
    };
    // A cohort's first row of each variant, and the spread of BETA and SE over its rows
    let cohort_agg_exprs = || -> Vec<Expr> {
        vec![
            col("BETA").first(),
            col("SE").first(),
            col("A1FREQ").first(),
            col("N").first(),
            col("ALLELE0").first(),
            col("ALLELE1").first(),
            col("BETA").min().alias("_beta_min"),
            col("BETA").max().alias("_beta_max"),
            col("SE").min().alias("_se_min"),
            col("SE").max().alias("_se_max"),
        ]
    };
    // A cohort's IVW terms, from its one row per variant
    let cohort_sum_exprs = || -> Vec<Expr> {
        vec![
            (lit(1.0f64) / (col("SE") * col("SE"))).alias("_w_sum"),
            (col("BETA") / (col("SE") * col("SE"))).alias("_wb_sum"),
            (col("BETA") * col("BETA") / (col("SE") * col("SE"))).alias("_wb2_sum"),
            (col("A1FREQ") * col("N")).alias("_afn_sum"),
            col("N").alias("_n_sum"),
            lit(1.0f64).alias("N_STUDIES"),
        ]
    };

    // Sums added to the running total
    let merge_agg_exprs = || -> Vec<Expr> {
        vec![
            col("_w_sum").sum(),
            col("_wb_sum").sum(),
            col("_wb2_sum").sum(),
            col("_afn_sum").sum(),
            col("_n_sum").sum(),
            col("N_STUDIES").sum(),
            col("ALLELE0").first(),
            col("ALLELE1").first(),
        ]
    };

    // ── 1. Accumulate the cohorts ────────────────────────────────────────────
    let t_accum = std::time::Instant::now();
    let mut accumulator: Option<DataFrame> = None;

    let polars_err = |e: PolarsError| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{e}"))
    };

    // Direction column of every cohort, and of those with files
    let mut dir_col_names: Vec<PlSmallStr> = Vec::new();
    let mut processed_dir_cols: Vec<PlSmallStr> = Vec::new();
    let cohort_label_str: String = cohort_files.iter().map(|(l, _)| l.as_str()).collect::<Vec<_>>().join(", ");

    for (cohort_idx, (label, files)) in cohort_files.iter().enumerate() {
        let dir_col_name_str = format!("_dir_{}", cohort_idx);
        let dir_col_name: PlSmallStr = dir_col_name_str.as_str().into();
        dir_col_names.push(dir_col_name.clone());

        if files.is_empty() { continue; }
        let t_c = std::time::Instant::now();
        processed_dir_cols.push(dir_col_name);

        let file_lfs: Vec<LazyFrame> = py.allow_threads(|| {
            files.iter()
                .map(|path| scan_regenie_df(path, sep, &cfg, &dtypes))
                .collect::<PolarsResult<Vec<_>>>()
        }).map_err(polars_err)?;

        let cohort_lf = if file_lfs.len() == 1 {
            file_lfs.into_iter().next().unwrap()
        } else {
            concat(
                file_lfs,
                UnionArgs { parallel: true, rechunk: false, to_supertypes: true, ..Default::default() },
            ).map_err(polars_err)?
        };

        // Copy each row into every bin its own A1FREQ falls in. Rows above the widest
        // bin are dropped first and the rest collected once, rather than rescanned per bin
        let cohort_lf = if bins.is_empty() {
            cohort_lf
        } else {
            let max_af = bins.iter().map(|(_, t)| *t).fold(0.0f64, f64::max);
            let rare = py.allow_threads(|| {
                cohort_lf.filter(col("A1FREQ").lt_eq(lit(max_af))).collect()
            }).map_err(polars_err)?;
            let per_bin: Vec<LazyFrame> = bins.iter().map(|(label, t)| {
                rare.clone().lazy()
                    .filter(col("A1FREQ").lt_eq(lit(*t)))
                    .with_column(lit(label.as_str()).alias("AAF_BIN"))
            }).collect();
            concat(per_bin, UnionArgs { parallel: true, rechunk: false, ..Default::default() })
                .map_err(polars_err)?
        };

        // One row per variant: a cohort's rows of a variant, from overlapping files or
        // repeated lines, count once if their BETA and SE agree
        let mut cohort_agg = cohort_agg_exprs();
        cohort_agg.push(col(dir_col_name_str.as_str()).first());
        let cohort_partial = cohort_lf
            .with_column(
                when(col("BETA").gt(lit(0.0f64)))
                    .then(lit("+"))
                    .otherwise(lit("-"))
                    .alias(dir_col_name_str.as_str()),
            )
            .group_by(group_cols())
            .agg(cohort_agg);

        let cohort_df = py.allow_threads(|| cohort_partial.collect())
            .map_err(polars_err)?;
        check_duplicates_agree(&cohort_df, label)?;
        let mut keep = group_cols();
        keep.extend([col("ALLELE0"), col("ALLELE1"), col(dir_col_name_str.as_str())]);
        keep.extend(cohort_sum_exprs());
        let cohort_df = py.allow_threads(|| cohort_df.lazy().select(keep).collect())
            .map_err(polars_err)?;

        accumulator = Some(match accumulator.take() {
            None => cohort_df,
            Some(prev) => {
                let merged = concat(
                    vec![prev.lazy(), cohort_df.lazy()],
                    // Diagonal, so other cohorts' direction columns are null
                    UnionArgs { diagonal: true, parallel: false, rechunk: false, to_supertypes: true, ..Default::default() },
                ).map_err(polars_err)?;
                let merge_group_cols = group_cols();
                let merge_dir_cols = processed_dir_cols.clone();
                py.allow_threads(|| {
                    let mut agg = merge_agg_exprs();
                    for dc in &merge_dir_cols {
                        // max() ignores nulls, so picks the one direction per variant
                        agg.push(col(dc.as_str()).max());
                    }
                    merged.group_by(merge_group_cols).agg(agg).collect()
                }).map_err(polars_err)?
            }
        });

        eprintln!(
            "[meta.rs/streaming] cohort {}/{} '{}': {:.3}s, accum rows: {}",
            cohort_idx + 1, n_cohorts, label,
            t_c.elapsed().as_secs_f64(),
            accumulator.as_ref().map_or(0, |df| df.height()),
        );
    }
    let accum_elapsed = t_accum.elapsed();

    let accum_df = accumulator.ok_or_else(|| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>("No data after accumulation")
    })?;

    // ── 2. Derive the statistics ─────────────────────────────────────────────
    let t_derive = std::time::Instant::now();
    let result_lf = accum_df.lazy()

        .with_columns([
            (col("_wb_sum") / col("_w_sum")).alias("BETA"),
            (lit(1.0f64) / col("_w_sum")).sqrt().alias("SE"),
            (col("_afn_sum") / col("_n_sum")).alias("A1FREQ"),
            col("_n_sum").cast(DataType::Int64).alias("N"),
            col("N_STUDIES").cast(DataType::UInt32),
        ])
        .with_columns([
            (col("BETA") / col("SE")).alias("Z"),
            when(
                (col("_wb2_sum") - col("_wb_sum") * col("_wb_sum") / col("_w_sum"))
                    .lt(lit(0.0f64)),
            )
            .then(lit(0.0f64))
            .otherwise(col("_wb2_sum") - col("_wb_sum") * col("_wb_sum") / col("_w_sum"))
            .alias("Q"),
        ])
        .with_columns([log10p_map_expr(), q_p_map_expr(), i2_map_expr()])
        .with_column({
            // '+' or '-' per cohort; '?' where the variant is missing or the cohort had no files
            let dir_exprs: Vec<Expr> = dir_col_names.iter().map(|dc| {
                if processed_dir_cols.contains(dc) {
                    col(dc.as_str()).fill_null(lit("?"))
                } else {
                    lit("?")
                }
            }).collect();
            concat_str(dir_exprs, "", false).alias("DIRECTIONS")
        });

    // ── 3. Select and sort ───────────────────────────────────────────────────
    let mut output_cols: Vec<Expr> = vec![
        col("CHROM"),
        col("GENPOS").alias("POS"),
        col("ID"),
        col("ALLELE0").alias("A0"),
        col("ALLELE1").alias("A1"),
        col("A1FREQ"),
        col("N"),
        col("TEST"),
        col("BETA"),
        col("SE"),
        col("Z"),
        col("LOG10P"),
        col("Q"),
        col("Q_P"),
        col("I2"),
        col("N_STUDIES"),
        col("DIRECTIONS"),
    ];
    let mut sort_cols = vec!["_chr_rank", "POS", "ID", "TEST"];
    if !bins.is_empty() {
        output_cols.push(col("AAF_BIN"));
        sort_cols.push("AAF_BIN");
    }

    let chr_rank = col("CHROM").map(
        |s| {
            let ca = s.as_materialized_series().str()?;
            let out: UInt32Chunked = ca.into_iter().map(|opt| {
                opt.map(|v| match v {
                    "X" => 23u32,
                    "Y" => 24,
                    "MT" => 25,
                    s => s.parse::<u32>().unwrap_or(26),
                })
            }).collect();
            Ok(out.into_column())
        },
        |_schema, field| Ok(Field::new(field.name().clone(), DataType::UInt32)),
    ).alias("_chr_rank");

    let mut result = py.allow_threads(|| {
        result_lf
            .with_columns([chr_rank])
            .select(
                output_cols.into_iter()
                    .chain(std::iter::once(col("_chr_rank")))
                    .collect::<Vec<_>>()
            )
            .sort(
                sort_cols,
                SortMultipleOptions::default().with_maintain_order(true),
            )
            .collect()
    }).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("compute error: {e}"))
    })?;

    result = result.drop("_chr_rank").map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("drop error: {e}"))
    })?;
    let derive_elapsed = t_derive.elapsed();
    let n_variants = result.height();

    // ── 4. Write ─────────────────────────────────────────────────────────────
    let t_write = std::time::Instant::now();
    let mut file = File::create(output_path)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{e}")))?;
    let write_err = |e: PolarsError| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("write error: {e}"));
    let io_err = |e: std::io::Error| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{e}"));

    match compression {
        "gz" => {
            let level = compression_level.map(|l| l as u32).unwrap_or(6);
            let mut enc = GzEncoder::new(file, GzCompression::new(level));
            writeln!(enc, "## Cohorts: {}", cohort_label_str).map_err(io_err)?;
            CsvWriter::new(enc)
                .with_separator(b'\t')
                .with_null_value("NA".into())
                .finish(&mut result)
                .map_err(write_err)?;
        }
        "zst" => {
            let level = compression_level.unwrap_or(3);
            let mut enc = zstd::Encoder::new(file, level)
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("zstd init: {e}")))?;
            writeln!(enc, "## Cohorts: {}", cohort_label_str).map_err(io_err)?;
            CsvWriter::new(&mut enc)
                .with_separator(b'\t')
                .with_null_value("NA".into())
                .finish(&mut result)
                .map_err(write_err)?;
            enc.finish()
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("zstd finish: {e}")))?;
        }
        _ => {
            writeln!(file, "## Cohorts: {}", cohort_label_str).map_err(io_err)?;
            CsvWriter::new(file)
                .with_separator(b'\t')
                .with_null_value("NA".into())
                .finish(&mut result)
                .map_err(write_err)?;
        }
    }
    let write_elapsed = t_write.elapsed();

    let n_files: usize = cohort_files.iter().map(|(_, f)| f.len()).sum();
    eprintln!(
        "[meta.rs/streaming] {}: accum={:.3}s  derive+sort={:.3}s  write={:.3}s  total={:.3}s  variants={}  cohorts={}  files={}",
        res_type,
        accum_elapsed.as_secs_f64(),
        derive_elapsed.as_secs_f64(),
        write_elapsed.as_secs_f64(),
        t_total.elapsed().as_secs_f64(),
        n_variants,
        n_cohorts,
        n_files,
    );

    Ok((n_variants, n_cohorts))
}
