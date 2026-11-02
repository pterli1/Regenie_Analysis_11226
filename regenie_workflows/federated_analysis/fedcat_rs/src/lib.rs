//! Rust extension of fedcat: meta-analysis, ACAT and annotation conversion.
//!
//! The Python API is declared, with its documentation, in `fedcat_rs.pyi`.

mod acat;
mod acat_core;
mod meta;
mod meta_merge;
mod streaming;

use pyo3::prelude::*;

/// Python module `fedcat_rs`.
#[pymodule]
fn fedcat_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(acat::acat_log10p_rs, m)?)?;
    m.add_function(wrap_pyfunction!(meta::z_to_log10p_rs, m)?)?;
    m.add_function(wrap_pyfunction!(meta::meta_analyse_regenie_streaming_rs, m)?)?;
    m.add_function(wrap_pyfunction!(meta_merge::meta_analyse_regenie_merge_rs, m)?)?;
    m.add_function(wrap_pyfunction!(streaming::compute_acat_streaming_rs, m)?)?;
    m.add_function(wrap_pyfunction!(streaming::convert_tsv_to_parquet_rs, m)?)?;
    Ok(())
}
