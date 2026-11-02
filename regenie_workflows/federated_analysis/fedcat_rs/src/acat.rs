//! Python function `acat_log10p_rs`: the ACAT of an array of -log10(P), computed by
//! [`acat_core::acat_log10p_core`].

use numpy::PyReadonlyArray1;
use pyo3::prelude::*;

use crate::acat_core;

/// Combine -log10(P) values by ACAT.
///
/// # Arguments
/// * `log10p_values` - -log10(P) values.
/// * `weights` - Weights of the same length, normalised to sum to 1; `None` for equal weights.
///
/// # Returns
/// -log10(P) of the combination, or NaN for no values.
///
/// # Errors
/// `ValueError` if `weights` is not the length of `log10p_values`.
#[pyfunction]
#[pyo3(signature = (log10p_values, weights=None))]
pub fn acat_log10p_rs<'py>(
    _py: Python<'py>,
    log10p_values: PyReadonlyArray1<'py, f64>,
    weights: Option<PyReadonlyArray1<'py, f64>>,
) -> PyResult<f64> {
    let vals = log10p_values.as_slice()?;
    let n = vals.len();
    if n == 0 {
        return Ok(f64::NAN);
    }

    let wts: Vec<f64> = match weights {
        None => vec![1.0 / n as f64; n],
        Some(w) => {
            let w_slice = w.as_slice()?;
            if w_slice.len() != n {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "weights and log10p_values must have same shape",
                ));
            }
            let sum: f64 = w_slice.iter().sum();
            w_slice.iter().map(|&x| x / sum).collect()
        }
    };

    Ok(acat_core::acat_log10p_core(vals, &wts))
}
