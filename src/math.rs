//! Numeric helpers shared by the algorithm.

use std::ops::Index;

/// Table of `ln(i!)` for `i` in `0..=n`, indexed by C `int`s.
///
/// The longest-path search can visit infeasible nodes and evaluate
/// `fact[-1]`. The C code then reads the word just before its workspace
/// buffer -- in practice malloc's chunk header, a denormal of order 1e-318 --
/// which is indistinguishable from 0 once added to a path length. Negative
/// indices therefore yield 0.0, reproducing the C results without the
/// out-of-bounds read.
pub(crate) struct LogFactorials(Vec<f64>);

impl LogFactorials {
    /// Builds the table with the C code's recurrence, which halves the number
    /// of `ln` calls (`ln(2k)! = ln(2k-1)! + ln 2 + ln k! - ln(k-1)!`).
    /// The exact recurrence matters for bit-identical results.
    pub(crate) fn new(n: i32) -> Self {
        let n = n as usize;
        let mut fact = vec![0f64; n + 1];
        if n >= 2 {
            fact[2] = 2f64.ln();
        }
        for i in (3..=n).step_by(2) {
            fact[i] = fact[i - 1] + (i as f64).ln();
            let j = i + 1;
            if j <= n {
                fact[j] = fact[i] + fact[2] + fact[j / 2] - fact[j / 2 - 1];
            }
        }
        LogFactorials(fact)
    }
}

impl Index<i32> for LogFactorials {
    type Output = f64;
    #[inline(always)]
    fn index(&self, i: i32) -> &f64 {
        if i < 0 {
            &0.0
        } else {
            &self.0[i as usize]
        }
    }
}

/// `ln(total! / prod(parts[i]!))`, the log of a multinomial coefficient
/// (C: `f9xact`). Subtracts in order, as the C code does.
pub(crate) fn log_multinomial(total: i32, parts: &[i32], fact: &LogFactorials) -> f64 {
    parts.iter().fold(fact[total], |d, &x| d - fact[x])
}

/// C's `fmax2`: `y` if `x < y`, else `x` (so NaN `y` yields `x`).
#[inline(always)]
pub(crate) fn c_fmax(x: f64, y: f64) -> f64 {
    if x < y {
        y
    } else {
        x
    }
}

/// C's `fmin2`: `x` if `x < y`, else `y`.
#[inline(always)]
pub(crate) fn c_fmin(x: f64, y: f64) -> f64 {
    if x < y {
        x
    } else {
        y
    }
}

/// Regularised lower incomplete gamma function P(p, y), i.e. the CDF of a
/// Gamma(p, 1) variable at `y` (C: `gammds`, Algorithm AS 147). Returns 0 for
/// non-positive arguments or when the leading factor underflows.
pub(crate) fn gamma_cdf(y: f64, p: f64) -> f64 {
    if y <= 0. || p <= 0. {
        return 0.;
    }
    let f = (p * y.ln() - log_gamma(p + 1.) - y).exp();
    if f == 0. {
        return 0.;
    }

    // Series expansion, summed until terms drop below 1e-6 relative.
    let mut c = 1.;
    let mut g = 1.;
    let mut a = p;
    loop {
        a += 1.;
        c *= y / a;
        g += c;
        // Negated so that a NaN term stops the loop, as in C's do-while.
        #[allow(clippy::neg_cmp_op_on_partial_ord)]
        if !(c > 1e-6 * g) {
            break;
        }
    }
    g * f
}

/// `ln Γ(x)` for `x >= 0` (C: `alogam`, ACM Algorithm 291): Stirling's series,
/// after shifting `x` up to at least 7. Returns 0 for negative `x`.
pub(crate) fn log_gamma(x: f64) -> f64 {
    const A1: f64 = 0.918938533204673;
    const A2: f64 = 5.95238095238e-4;
    const A3: f64 = 7.93650793651e-4;
    const A4: f64 = 0.002777777777778;
    const A5: f64 = 0.083333333333333;

    if x < 0. {
        return 0.;
    }
    let mut y = x;
    let mut f = 0.;
    if y < 7. {
        // f = x (x+1) ... (y-1) with y the first value >= 7.
        f = y;
        loop {
            y += 1.;
            if y >= 7. {
                break;
            }
            f *= y;
        }
        f = -f.ln();
    }
    let z = 1. / (y * y);
    f + (y - 0.5) * y.ln() - y + A1 + (((-A2 * z + A3) * z - A4) * z + A5) / y
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_factorials() {
        let fact = LogFactorials::new(20);
        let mut exact = 0f64;
        for i in 1..=20 {
            exact += (i as f64).ln();
            assert!((fact[i] - exact).abs() < 1e-12, "{i}");
        }
        assert_eq!(fact[-1], 0.0);
    }

    #[test]
    fn gamma() {
        assert!((log_gamma(5.) - 24f64.ln()).abs() < 1e-9); // ACM 291 is accurate to ~1e-10
                                                            // P(1, y) = 1 - e^-y
        assert!((gamma_cdf(0.7, 1.) - (1. - (-0.7f64).exp())).abs() < 1e-6);
        assert_eq!(gamma_cdf(0., 2.), 0.);
    }
}
