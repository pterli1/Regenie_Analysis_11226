//! ACAT (aggregated Cauchy association test) in -log10(P) space, used by `acat.rs` and
//! `streaming.rs`.

use std::f64::consts::{LN_10, PI};

/// -log10(P) from which a term is accumulated in log space.
pub const ASYMPTOTIC_LOG10P: f64 = 15.0;

/// ln(π)
const LN_PI: f64 = 1.144_729_885_849_400_2;
/// ln(T) above which the asymptotic p ≈ 1/(πT) is used, as in REGENIE (exp(16) ≈ 8.886 × 10⁶).
const LN_ASYMPTOTIC_T: f64 = 16.0;

/// Convert ln(T) of a positive ACAT statistic to -log10(P): exact Cauchy survival, or
/// p ≈ 1/(πT) for large T.
#[inline]
fn lsum_to_log10p(lsum: f64) -> f64 {
    if lsum >= LN_ASYMPTOTIC_T {
        (lsum + LN_PI) / LN_10
    } else {
        let t = lsum.exp();
        let p = 0.5 - t.atan() / PI;
        if p <= 0.0 { 308.0 } else { -p.log10() }
    }
}

/// Combine -log10(P) values by ACAT, as REGENIE's `get_acat_robust` does.
///
/// Terms with -log10(P) >= [`ASYMPTOTIC_LOG10P`] are summed in log space (T_A), so do
/// not overflow; the others in linear space (T_B, possibly negative). The two are
/// combined with `ln_1p`, whatever their signs.
///
/// # Arguments
/// * `vals` - -log10(P) values.
/// * `wts` - Weights of the same length, summing to 1.
///
/// # Returns
/// -log10(P) of the combination, at most 308 unless T_A is present, or NaN for no values.
#[inline]
pub fn acat_log10p_core(vals: &[f64], wts: &[f64]) -> f64 {
    let n = vals.len();
    if n == 0 {
        return f64::NAN;
    }

    let mut ln_ta: Vec<f64> = Vec::new();
    let mut tb: f64 = 0.0;

    for k in 0..n {
        if vals[k] >= ASYMPTOTIC_LOG10P {
            // ln(w_k × 10^vals[k] / π) = ln(w_k) + vals[k]·ln(10) − ln(π)
            ln_ta.push(
                wts[k].max(f64::MIN_POSITIVE).ln() + vals[k] * LN_10 - LN_PI,
            );
        } else {
            // p < 1, as p = 1 would give tan(−π/2) = −∞
            let p_k = f64::powf(10.0, -vals[k]).min(0.999);
            // tan((0.5 − p)π) as 1/tan(pπ): forming 0.5 − p first loses precision for small p
            tb += wts[k] / (p_k * PI).tan();
        }
    }

    let has_a = !ln_ta.is_empty();

    // --- Only T_B, as in REGENIE ---
    if !has_a {
        if tb >= 8_886_111.0_f64 {
            return (tb.ln() + LN_PI) / LN_10;
        }
        let p = 0.5 - tb.atan() / PI;
        if p <= 0.0 { return 308.0; }
        if p >= 1.0 { return 0.0; }
        return -p.log10();
    }

    // --- ln(T_A) by log-sum-exp ---
    let max_ln = ln_ta.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let l_ta = max_ln + ln_ta.iter().map(|&x| (x - max_ln).exp()).sum::<f64>().ln();

    // --- Only T_A ---
    if tb == 0.0 {
        return lsum_to_log10p(l_ta);
    }

    // --- Both ---
    if tb < 0.0 {
        let l_abs_tb = tb.abs().ln();
        if l_abs_tb < l_ta {
            // T_total = T_A − |T_B| > 0
            let lsum = l_ta + (-(l_abs_tb - l_ta).exp()).ln_1p();
            lsum_to_log10p(lsum)
        } else {
            // T_total ≤ 0 → p ≥ 0.5, log10p near 0
            let lsum = l_abs_tb + (-(l_ta - l_abs_tb).exp()).ln_1p();
            if lsum >= LN_ASYMPTOTIC_T {
                // p ≈ 1 − 1/(π|T|)
                let lp = (-(-(lsum + LN_PI)).exp()).ln_1p();
                (-lp / LN_10).max(0.0)
            } else {
                let t_neg = lsum.exp();
                let p = 0.5 + t_neg.atan() / PI;
                if p >= 1.0 { 0.0 } else { -p.log10() }
            }
        }
    } else {
        // T_B > 0: T_total = T_A + T_B > 0
        let l_tb = tb.ln();
        let lsum = l_ta.max(l_tb) + (-(l_ta - l_tb).abs()).exp().ln_1p();
        lsum_to_log10p(lsum)
    }
}

/// Unnormalised ACAT weight of a variant, as REGENIE and GCTA weight them.
///
/// The weight is maf × (1 − maf) × Beta(maf; 1, 25)², which is proportional to
/// maf × (1 − maf)^49 ([Liu et al. 2019](https://doi.org/10.1016/j.ajhg.2019.01.002)).
/// Weights must be normalised to sum to 1 within each combination.
///
/// # Arguments
/// * `af` - Allele frequency; the minor allele frequency is taken from it.
#[inline]
pub fn weight_standard(af: f64) -> f64 {
    let maf = if af <= 0.5 { af } else { 1.0 - af };
    maf * (1.0 - maf).powi(49)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_acat_single_small_p_is_exact() {
        // A single p-value combines to itself, including near ASYMPTOTIC_LOG10P
        for v in [12.0, 13.0, 14.0, 14.9] {
            let result = acat_log10p_core(&[v], &[1.0]);
            assert!((result - v).abs() < 1e-6, "{v}: {result}");
        }
    }

    #[test]
    fn test_acat_empty() {
        assert!(acat_log10p_core(&[], &[]).is_nan());
    }

    #[test]
    fn test_acat_single_value() {
        let result = acat_log10p_core(&[5.0], &[1.0]);
        assert!(result > 0.0 && result.is_finite());
    }

    #[test]
    fn test_acat_extreme_ta_only() {
        // A tiny weight on a huge -log10(P) must not overflow (10^350) in T_A
        let vals = vec![350.0, 1.0];
        let wts = vec![1e-100_f64 / (1e-100 + 1.0), 1.0 / (1e-100 + 1.0)];
        let result = acat_log10p_core(&vals, &wts);
        assert!(result.is_finite(), "overflow bug: got {result}");
    }

    #[test]
    fn test_acat_extreme_dominant() {
        // Very large logp dominates; result should be close to vals[0]
        let vals = vec![400.0, 1.0];
        let wts = vec![0.5, 0.5];
        let result = acat_log10p_core(&vals, &wts);
        assert!(result > 300.0 && result.is_finite());
    }

    #[test]
    fn test_acat_pval_one_no_nan() {
        // -log10(P) = 0 is p = 1, where tan(-π/2) = -∞
        let vals = vec![0.0, 3.0];
        let wts = vec![0.5, 0.5];
        let result = acat_log10p_core(&vals, &wts);
        assert!(result.is_finite() && result >= 0.0, "p=1 bug: got {result}");
    }

    #[test]
    fn test_acat_tb_negative_no_asymptotic() {
        // Only weak p-values: T_B < 0 and no T_A
        let vals = vec![0.1, 0.05, 0.02];
        let wts = vec![1.0 / 3.0; 3];
        let result = acat_log10p_core(&vals, &wts);
        assert!(result.is_finite() && result >= 0.0, "tb<0 no_a bug: got {result}");
        assert!(result < 1.0);
    }
}
