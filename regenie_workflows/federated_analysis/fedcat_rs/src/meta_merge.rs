//! IVW meta-analysis by a K-way merge of sorted inputs, in memory proportional to the
//! number of files.
//!
//! All files are read at once, merged on `(chromosome, GENPOS, ID, TEST)`, and each
//! variant written as soon as all its rows are read. Input must be sorted, as REGENIE
//! writes it; out-of-order rows are warned about.

use pyo3::prelude::*;
use pyo3::types::PyDict;
use statrs::distribution::{ChiSquared, ContinuousCDF};

use flate2::read::GzDecoder;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use crate::meta::{parse_aaf_bins, z_to_log10p_scalar, ColumnConfig, DUPLICATE_TOLERANCE};

// ── Row representation ───────────────────────────────────────────────────────

/// A valid row of a REGENIE file, with the fields the meta-analysis needs.
struct RegenieRow {
    chrom: String,
    chrom_rank: u32,
    genpos: i64,
    id: String,
    test: String,
    a0: String,
    a1: String,
    a1freq: f64,
    n: f64,
    beta: f64,
    se: f64,
}

/// Merge order of rows: chromosome rank, position, ID, test.
#[derive(Eq, PartialEq)]
struct RowKey {
    chrom_rank: u32,
    genpos: i64,
    id: String,
    test: String,
}

impl Ord for RowKey {
    fn cmp(&self, other: &Self) -> Ordering {
        self.chrom_rank.cmp(&other.chrom_rank)
            .then(self.genpos.cmp(&other.genpos))
            .then(self.id.cmp(&other.id))
            .then(self.test.cmp(&other.test))
    }
}

impl PartialOrd for RowKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl RegenieRow {
    /// The row's merge key.
    fn key(&self) -> RowKey {
        RowKey {
            chrom_rank: self.chrom_rank,
            genpos: self.genpos,
            id: self.id.clone(),
            test: self.test.clone(),
        }
    }
}

// ── Per-file reader ──────────────────────────────────────────────────────────

/// Positions of the needed columns in a file's header.
struct ColIndices {
    chrom: usize,
    genpos: usize,
    id: usize,
    a0: usize,
    a1: usize,
    a1freq: usize,
    n: usize,
    test: usize,
    beta: usize,
    se: usize,
    extra: Option<usize>,
}

/// Reader of the rows of one REGENIE file, plain, gzip or zstd.
struct FileReader {
    /// Decompressed line reader.
    reader: Box<dyn BufRead + Send>,
    /// Column positions.
    idx: ColIndices,
    /// Field separator.
    sep: u8,
    /// Columns in the header.
    n_cols: usize,
    /// Position of the file among all inputs.
    #[allow(dead_code)]
    file_id: usize,
    /// Key, chromosome and position of the last row returned, to check the sort order.
    last_key: Option<(RowKey, String, i64)>,
    /// File path, for messages.
    path: String,
    /// Line buffer, reused for every line.
    line_buf: String,
}

/// Sort rank of a chromosome: 1-22, X (or 23), Y, MT, then others.
fn chrom_rank(c: &str) -> u32 {
    match c {
        "X" | "23" => 23,
        "Y" => 24,
        "MT" => 25,
        s => s.parse::<u32>().unwrap_or(26),
    }
}

impl FileReader {
    /// Open a file and read its header, skipping `##` lines.
    ///
    /// # Errors
    /// `IOError` if it cannot be read; `ValueError` if it has no header or lacks a column.
    fn open(path: &str, sep: u8, cfg: &ColumnConfig, file_id: usize) -> PyResult<Self> {
        let file = File::open(path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{e}")))?;
        let mut reader: Box<dyn BufRead + Send> = if path.ends_with(".gz") {
            Box::new(BufReader::with_capacity(1 << 18, GzDecoder::new(file)))
        } else if path.ends_with(".zst") || path.ends_with(".zstd") {
            Box::new(BufReader::with_capacity(1 << 18,
                zstd::Decoder::new(file)
                    .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("zstd: {e}")))?))
        } else {
            Box::new(BufReader::with_capacity(1 << 18, file))
        };

        let mut header_line = String::new();
        loop {
            header_line.clear();
            let bytes = reader.read_line(&mut header_line)
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{e}")))?;
            if bytes == 0 {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                    format!("No header found in {path}"),
                ));
            }
            let trimmed = header_line.trim_start();
            if trimmed.starts_with("##") || trimmed.trim().is_empty() {
                continue;
            }
            break;
        }
        let headers: Vec<&str> = header_line.trim_end().split(sep as char).collect();
        let n_cols = headers.len();

        let find = |name: &str| -> PyResult<usize> {
            headers.iter().position(|h| *h == name)
                .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyValueError, _>(
                    format!("column '{}' not in header of {}", name, path)))
        };

        let idx = ColIndices {
            chrom: find(&cfg.chrom)?,
            genpos: find(&cfg.genpos)?,
            id: find(&cfg.id)?,
            a0: find(&cfg.allele0)?,
            a1: find(&cfg.allele1)?,
            a1freq: find(&cfg.a1freq)?,
            n: find(&cfg.n)?,
            test: find(&cfg.test)?,
            beta: find(&cfg.beta)?,
            se: find(&cfg.se)?,
            extra: headers.iter().position(|h| *h == cfg.extra.as_deref().unwrap_or("EXTRA")),
        };

        Ok(Self { reader, idx, sep, n_cols, file_id, last_key: None, path: path.to_string(), line_buf: String::with_capacity(256) })
    }

    /// Read the next valid row, or `None` at the end of the file.
    ///
    /// Skips short rows, `TEST_FAIL` rows, and rows where SE <= 0 or BETA, SE, A1FREQ
    /// or N is not finite; chromosome 23 becomes X.
    fn next_row(&mut self) -> PyResult<Option<RegenieRow>> {
        let sep_char = self.sep as char;
        loop {
            self.line_buf.clear();
            let bytes = self.reader.read_line(&mut self.line_buf)
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{e}")))?;
            if bytes == 0 { return Ok(None); }

            let trimmed = self.line_buf.trim_end();
            if trimmed.is_empty() || trimmed.starts_with("##") {
                continue;
            }
            let fields: Vec<&str> = trimmed.splitn(self.n_cols + 1, sep_char).collect();
            if fields.len() < self.n_cols { continue; }

            if let Some(ei) = self.idx.extra {
                if fields.get(ei).map_or(false, |e| *e == "TEST_FAIL") {
                    continue;
                }
            }

            // Parsing accepts "nan" and "inf", so finiteness is checked too
            let se: f64 = match fields[self.idx.se].parse() {
                Ok(v) if v > 0.0 && f64::is_finite(v) => v,
                _ => continue,
            };
            let beta: f64 = match fields[self.idx.beta].parse() {
                Ok(v) if f64::is_finite(v) => v,
                _ => continue,
            };
            let a1freq: f64 = match fields[self.idx.a1freq].parse() {
                Ok(v) if f64::is_finite(v) => v,
                _ => continue,
            };
            let n: f64 = match fields[self.idx.n].parse() {
                Ok(v) if f64::is_finite(v) => v,
                _ => continue,
            };

            let chrom_raw = fields[self.idx.chrom];
            let chrom = if chrom_raw == "23" { "X".to_string() } else { chrom_raw.to_string() };
            let cr = chrom_rank(&chrom);
            let genpos: i64 = fields[self.idx.genpos].parse().unwrap_or(0);
            let id = fields[self.idx.id].to_string();
            let test = fields[self.idx.test].to_string();

            let row = RegenieRow {
                chrom,
                chrom_rank: cr,
                genpos,
                id,
                test,
                a0: fields[self.idx.a0].to_string(),
                a1: fields[self.idx.a1].to_string(),
                a1freq,
                n,
                beta,
                se,
            };

            // Out-of-order rows are warned about, not rejected
            let key = row.key();
            if let Some((ref last, ref last_chrom, last_pos)) = self.last_key {
                if key < *last {
                    eprintln!(
                        "[meta_merge] WARNING: file {} not sorted at pos {}:{} (prev {}:{})",
                        self.path, row.chrom, row.genpos, last_chrom, last_pos
                    );
                }
            }
            self.last_key = Some((key, row.chrom.clone(), row.genpos));

            return Ok(Some(row));
        }
    }
}

// ── Min-heap entry for k-way merge ───────────────────────────────────────────

/// A file's next row in the merge heap, ordered in reverse so the max-heap pops the lowest key.
struct HeapEntry {
    /// The row.
    row: RegenieRow,
    /// File it came from.
    file_id: usize,
}

impl Eq for HeapEntry {}
impl PartialEq for HeapEntry {
    fn eq(&self, other: &Self) -> bool {
        self.row.key() == other.row.key() && self.file_id == other.file_id
    }
}
impl Ord for HeapEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        other.row.key().cmp(&self.row.key())
            .then(other.file_id.cmp(&self.file_id))
    }
}
impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

// ── IVW accumulator ──────────────────────────────────────────────────────────

/// IVW sums of one variant (and bin) across cohorts.
struct IvwAccum {
    chrom: String,
    genpos: i64,
    id: String,
    test: String,
    /// Reference allele of the first cohort.
    a0: String,
    /// Tested allele of the first cohort.
    a1: String,
    /// Sum of 1/SE².
    w_sum: f64,
    /// Sum of BETA/SE².
    wb_sum: f64,
    /// Sum of BETA²/SE².
    wb2_sum: f64,
    /// Sum of A1FREQ × N.
    afn_sum: f64,
    /// Sum of N.
    n_sum: f64,
    /// Cohorts added.
    n_studies: u32,
    /// `(cohort, beta, se)` of each cohort added, to spot one cohort's duplicate rows.
    cohort_contribs: Vec<(usize, f64, f64)>,
}

impl IvwAccum {
    /// Start the sums with a cohort's row.
    fn new(row: &RegenieRow, cohort_id: usize) -> Self {
        let w = 1.0 / (row.se * row.se);
        Self {
            chrom: row.chrom.clone(),
            genpos: row.genpos,
            id: row.id.clone(),
            test: row.test.clone(),
            a0: row.a0.clone(),
            a1: row.a1.clone(),
            w_sum: w,
            wb_sum: w * row.beta,
            wb2_sum: w * row.beta * row.beta,
            afn_sum: row.a1freq * row.n,
            n_sum: row.n,
            n_studies: 1,
            cohort_contribs: vec![(cohort_id, row.beta, row.se)],
        }
    }

    /// Add a cohort's row.
    fn add(&mut self, row: &RegenieRow, cohort_id: usize) {
        let w = 1.0 / (row.se * row.se);
        self.w_sum += w;
        self.wb_sum += w * row.beta;
        self.wb2_sum += w * row.beta * row.beta;
        self.afn_sum += row.a1freq * row.n;
        self.n_sum += row.n;
        self.n_studies += 1;
        self.cohort_contribs.push((cohort_id, row.beta, row.se));
    }

    /// The variant's merge key.
    #[allow(dead_code)]
    fn key(&self) -> RowKey {
        RowKey {
            chrom_rank: chrom_rank(&self.chrom),
            genpos: self.genpos,
            id: self.id.clone(),
            test: self.test.clone(),
        }
    }
}

/// P-value of Cochran's Q on `ns - 1` degrees of freedom; `None` (written NA) for one study.
fn cochran_q_pvalue(q: f64, ns: u32) -> Option<f64> {
    if ns <= 1 { return None; }
    let df = (ns - 1) as f64;
    Some(ChiSquared::new(df).unwrap().sf(q))
}

/// I² heterogeneity index, (Q - df) / Q clipped to [0, 1]; 0 for one study.
#[inline]
fn i2_stat(q: f64, ns: u32) -> f64 {
    if ns <= 1 { return 0.0; }
    let df = (ns - 1) as f64;
    if q > 0.0 { ((q - df) / q).clamp(0.0, 1.0) } else { 0.0 }
}

/// Write a variant's output row from its sums, ending with the `AAF_BIN` column when binning.
fn flush_accum(acc: &IvwAccum, w: &mut impl Write, n_cohorts: usize, bin: Option<&str>) -> std::io::Result<()> {
    let beta = acc.wb_sum / acc.w_sum;
    let se = (1.0 / acc.w_sum).sqrt();
    let z = beta / se;
    let log10p = z_to_log10p_scalar(z);
    let a1f = acc.afn_sum / acc.n_sum;
    let n_i = acc.n_sum as i64;
    let q_raw = acc.wb2_sum - acc.wb_sum * acc.wb_sum / acc.w_sum;
    let q = if q_raw < 0.0 { 0.0 } else { q_raw };
    let q_p = cochran_q_pvalue(q, acc.n_studies).map_or_else(|| "NA".to_string(), |p| format!("{p:.4e}"));
    let i2 = i2_stat(q, acc.n_studies);

    let mut dirs = vec![b'?'; n_cohorts];
    for &(cohort_id, beta_c, _) in &acc.cohort_contribs {
        if cohort_id < n_cohorts {
            dirs[cohort_id] = if beta_c > 0.0 { b'+' } else { b'-' };
        }
    }
    // Safety: dirs holds only the ASCII bytes '+', '-' and '?'
    let dirs_str = unsafe { String::from_utf8_unchecked(dirs) };

    write!(
        w,
        "{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.4}\t{:.4}\t{}\t{:.4}\t{}\t{}",
        acc.chrom, acc.genpos, acc.id, acc.a0, acc.a1,
        a1f, n_i, acc.test,
        beta, se, z, log10p, q, q_p, i2, acc.n_studies, dirs_str
    )?;
    match bin {
        Some(label) => writeln!(w, "\t{label}"),
        None => writeln!(w),
    }
}

// ── Public entry point ───────────────────────────────────────────────────────

/// Meta-analyse sorted REGENIE results by IVW fixed effects, merging all files at once.
///
/// Output is as for [`crate::meta::meta_analyse_regenie_streaming_rs`], including a
/// cohort's rows of one variant counting once if their BETA and SE agree.
///
/// # Arguments
/// As for [`crate::meta::meta_analyse_regenie_streaming_rs`]; input files must be
/// sorted by chromosome, GENPOS, ID and TEST.
///
/// # Returns
/// `(n_rows_written, n_cohorts)`.
///
/// # Errors
/// `ValueError` for no input files, a missing column, or a cohort giving one variant
/// different results; `IOError` if a file cannot be read or written.
#[pyfunction]
#[pyo3(signature = (cohort_files, output_path, res_type, column_config=None, compression="zst", compression_level=None, aaf_bins=None))]
pub fn meta_analyse_regenie_merge_rs(
    _py: Python<'_>,
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

    // Every file, with its cohort
    let mut all_files: Vec<String> = Vec::new();
    let mut file_to_cohort: Vec<usize> = Vec::new();
    for (cohort_id, (_label, files)) in cohort_files.iter().enumerate() {
        for f in files {
            file_to_cohort.push(cohort_id);
            all_files.push(f.clone());
        }
    }

    if all_files.is_empty() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("No input files"));
    }

    let t_total = std::time::Instant::now();

    // Configured separator, or tab if the first file's header has one, else space
    let sep = if let Some(s) = cfg.separator { s as u8 } else {
        let probe_file = File::open(&all_files[0])
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{e}")))?;
        let probe_reader: Box<dyn BufRead> = if all_files[0].ends_with(".gz") {
            Box::new(BufReader::new(GzDecoder::new(probe_file)))
        } else {
            Box::new(BufReader::new(probe_file))
        };
        let mut first_line = String::new();
        let mut br = BufReader::new(probe_reader);
        loop {
            first_line.clear();
            let bytes = br.read_line(&mut first_line).unwrap_or(0);
            if bytes == 0 { break; }
            let trimmed = first_line.trim_start();
            if trimmed.starts_with("##") || trimmed.trim().is_empty() {
                continue;
            }
            break;
        }
        if first_line.contains('\t') { b'\t' } else { b' ' }
    };

    // Open every file, with its first row in the heap
    let k = all_files.len();
    let mut readers: Vec<FileReader> = Vec::with_capacity(k);
    let mut heap: BinaryHeap<HeapEntry> = BinaryHeap::with_capacity(k);

    for (i, path) in all_files.iter().enumerate() {
        let mut fr = FileReader::open(path, sep, &cfg, i)?;
        if let Some(row) = fr.next_row()? {
            heap.push(HeapEntry { row, file_id: i });
        }
        readers.push(fr);
    }

    let file = File::create(output_path)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{e}")))?;

    let writer: Box<dyn Write> = match compression {
        "gz" => {
            let level = compression_level.map(|l| l as u32).unwrap_or(6);
            Box::new(flate2::write::GzEncoder::new(file, flate2::Compression::new(level)))
        }
        "zst" => {
            let level = compression_level.unwrap_or(3);
            let enc = zstd::Encoder::new(file, level)
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("zstd: {e}")))?;
            Box::new(enc.auto_finish())
        }
        _ => Box::new(file),
    };
    let mut bw = BufWriter::with_capacity(1 << 20, writer);

    let cohort_names: Vec<&str> = cohort_files.iter().map(|(l, _)| l.as_str()).collect();
    writeln!(bw, "## Cohorts: {}", cohort_names.join(", "))
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{e}")))?;
    let bin_header = if bins.is_empty() { "" } else { "\tAAF_BIN" };
    writeln!(bw, "CHROM\tPOS\tID\tA0\tA1\tA1FREQ\tN\tTEST\tBETA\tSE\tZ\tLOG10P\tQ\tQ_P\tI2\tN_STUDIES\tDIRECTIONS{bin_header}")
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{e}")))?;

    // One accumulator per bin, or one without bins; a row joins every bin its A1FREQ falls in
    let n_slots = bins.len().max(1);
    let in_slot = |slot: usize, a1freq: f64| bins.is_empty() || a1freq <= bins[slot].1;
    let slot_label = |slot: usize| bins.get(slot).map(|(l, _)| l.as_str());

    let mut n_variants: usize = 0;
    let mut current: Vec<Option<IvwAccum>> = (0..n_slots).map(|_| None).collect();
    let mut current_key: Option<(String, i64, String, String)> = None;
    let io_err = |e: std::io::Error| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{e}"));

    // K-way merge loop
    while let Some(entry) = heap.pop() {
        let HeapEntry { row, file_id } = entry;
        let cohort_id = file_to_cohort[file_id];

        // Replace the row with the next from its file
        if let Some(next_row) = readers[file_id].next_row()? {
            heap.push(HeapEntry { row: next_row, file_id });
        }

        let same_key = current_key.as_ref().map_or(false, |(c, g, i, t)| {
            *i == row.id && *t == row.test && *c == row.chrom && *g == row.genpos
        });
        if !same_key {
            for (slot, acc) in current.iter_mut().enumerate() {
                if let Some(acc) = acc.take() {
                    flush_accum(&acc, &mut bw, n_cohorts, slot_label(slot)).map_err(io_err)?;
                    n_variants += 1;
                }
            }
            current_key = Some((row.chrom.clone(), row.genpos, row.id.clone(), row.test.clone()));
        }

        for slot in 0..n_slots {
            if !in_slot(slot, row.a1freq) {
                continue;
            }
            match &mut current[slot] {
                Some(acc) => {
                    // A cohort's second row of a variant must match the first
                    if let Some((_, prev_beta, prev_se)) = acc.cohort_contribs.iter().find(|(c, _, _)| *c == cohort_id) {
                        let beta_match = (prev_beta - row.beta).abs() <= DUPLICATE_TOLERANCE;
                        let se_match = (prev_se - row.se).abs() <= DUPLICATE_TOLERANCE;
                        if beta_match && se_match {
                            continue;
                        }
                        let cohort_label = &cohort_files[cohort_id].0;
                        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                            format!(
                                "Conflicting results for variant '{}' (test={}) in cohort '{}': \
                                 BETA={:.6}/SE={:.6} vs earlier BETA={:.6}/SE={:.6}",
                                row.id, row.test, cohort_label,
                                row.beta, row.se, prev_beta, prev_se,
                            )
                        ));
                    }
                    acc.add(&row, cohort_id);
                }
                None => current[slot] = Some(IvwAccum::new(&row, cohort_id)),
            }
        }
    }

    for (slot, acc) in current.iter().enumerate() {
        if let Some(acc) = acc {
            flush_accum(acc, &mut bw, n_cohorts, slot_label(slot)).map_err(io_err)?;
            n_variants += 1;
        }
    }

    bw.flush().map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("{e}")))?;

    eprintln!(
        "[meta.rs/merge] {}: total={:.3}s  variants={}  cohorts={}  files={}",
        res_type, t_total.elapsed().as_secs_f64(), n_variants, n_cohorts, k,
    );

    Ok((n_variants, n_cohorts))
}
