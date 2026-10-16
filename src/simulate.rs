//! Monte Carlo p-values (as R's `fisher.test(simulate.p.value = TRUE)`).
//!
//! Draws `samples` random tables with the observed margins and counts those
//! at most as probable as the observed table: `p = (1 + hits) / (samples +
//! 1)`. Tables are drawn cell by cell: given what is left of the row and of
//! the columns, a cell is hypergeometric (Patefield, AS 159), sampled by
//! inversion outward from its mode.
//!
//! Samples are split into fixed chunks, each with its own random stream
//! derived from the seed, so the result depends on the seed but not on the
//! number of threads.

use std::hash::{BuildHasher, RandomState};
use std::num::NonZeroUsize;
use std::sync::atomic::{AtomicU64, Ordering};

use crate::math::LogFactorials;
use crate::{Error, Table};

/// Samples per random stream.
const CHUNK: u64 = 256;

/// R's tolerance for ties with the observed table (`1 + 64 * eps`).
const ALMOST_ONE: f64 = 1. + 64. * f64::EPSILON;

/// Settings for [`fisher_simulated`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Simulation {
    /// Number of random tables. Default 10 000.
    pub samples: u64,
    /// Seed of the random streams; `None` picks a random one. Default
    /// `None`.
    pub seed: Option<u64>,
    /// Number of threads; `None` uses all available cores. Default `None`.
    pub threads: Option<NonZeroUsize>,
}

impl Default for Simulation {
    fn default() -> Self {
        Simulation {
            samples: 10_000,
            seed: None,
            threads: None,
        }
    }
}

/// Result of [`fisher_simulated`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Simulated {
    /// `(1 + at_least_as_extreme) / (samples + 1)`.
    pub p_value: f64,
    pub samples: u64,
    /// Random tables at most as probable as the observed one.
    pub at_least_as_extreme: u64,
}

impl Simulated {
    /// Estimated standard error of `p_value`.
    pub fn standard_error(&self) -> f64 {
        let p = self.p_value;
        (p * (1. - p) / (self.samples + 1) as f64).sqrt()
    }
}

/// Monte Carlo estimate of Fisher's exact test p-value for `table`.
///
/// Unlike [`fisher_exact`](crate::fisher_exact) this has no workspace
/// limit: time grows with `samples` times the number of cells (and, for
/// huge counts, with their square root).
pub fn fisher_simulated(table: &Table, simulation: &Simulation) -> Result<Simulated, Error> {
    let (nrow, ncol) = (table.nrows(), table.ncols());
    let mut ntot: i64 = 0;
    for i in 0..nrow {
        for j in 0..ncol {
            let x = table.get(i, j);
            if x < 0 {
                return Err(Error::NegativeEntry);
            }
            ntot = ntot.checked_add(x).ok_or(Error::IntegerOverflow)?;
        }
    }
    if ntot == 0 {
        return Err(Error::AllZero);
    }

    // Empty rows and columns are always empty: drop them.
    let rows: Vec<i64> = (0..nrow)
        .map(|i| (0..ncol).map(|j| table.get(i, j)).sum())
        .filter(|&r| r > 0)
        .collect();
    let cols: Vec<i64> = (0..ncol)
        .map(|j| (0..nrow).map(|i| table.get(i, j)).sum())
        .filter(|&c| c > 0)
        .collect();
    let samples = simulation.samples;
    if rows.len() < 2 || cols.len() < 2 {
        // Only one table has these margins.
        return Ok(Simulated {
            p_value: 1.,
            samples,
            at_least_as_extreme: samples,
        });
    }

    let fact = LnFactorial::new(ntot);
    let observed = (0..nrow)
        .flat_map(|i| (0..ncol).map(move |j| (i, j)))
        .fold(0., |s, (i, j)| s + fact.get(table.get(i, j)));
    let sampler = Sampler {
        rows: &rows,
        cols: &cols,
        ntot,
        fact: &fact,
        // A table is at most as probable as the observed one when its
        // `sum ln(cell!)` is at least as large (with R's tolerance).
        threshold: observed / ALMOST_ONE,
    };

    let seed = simulation
        .seed
        .unwrap_or_else(|| RandomState::new().hash_one(0u64));
    let chunks = samples.div_ceil(CHUNK);
    let threads = simulation
        .threads
        .or_else(|| std::thread::available_parallelism().ok())
        .map_or(1, NonZeroUsize::get)
        .min(usize::try_from(chunks).unwrap_or(usize::MAX));

    let next = AtomicU64::new(0);
    let work = || {
        let mut col_left = vec![0; cols.len()];
        let mut hits = 0;
        loop {
            let chunk = next.fetch_add(1, Ordering::Relaxed);
            if chunk >= chunks {
                return hits;
            }
            let n = CHUNK.min(samples - chunk * CHUNK);
            let mut rng = Xoshiro256pp::new(seed, chunk);
            hits += (0..n)
                .filter(|_| sampler.extreme(&mut rng, &mut col_left))
                .count() as u64;
        }
    };
    let hits = if threads <= 1 {
        work()
    } else {
        std::thread::scope(|s| {
            let workers: Vec<_> = (0..threads).map(|_| s.spawn(work)).collect();
            workers
                .into_iter()
                .map(|w| w.join().unwrap_or_else(|p| std::panic::resume_unwind(p)))
                .sum()
        })
    };

    Ok(Simulated {
        p_value: (hits as f64 + 1.) / (samples as f64 + 1.),
        samples,
        at_least_as_extreme: hits,
    })
}

/// Draws random tables with fixed margins.
struct Sampler<'a> {
    rows: &'a [i64],
    cols: &'a [i64],
    ntot: i64,
    fact: &'a LnFactorial,
    threshold: f64,
}

impl Sampler<'_> {
    /// Draws a table and tells whether it is at least as extreme as the
    /// observed one. `col_left` is scratch space of `cols.len()`.
    fn extreme(&self, rng: &mut Xoshiro256pp, col_left: &mut [i64]) -> bool {
        let fact = self.fact;
        let rows = &self.rows[..self.rows.len() - 1];
        col_left.copy_from_slice(self.cols);
        let mut total_left = self.ntot;
        let mut stat = 0.;
        for &r in rows {
            // Spread the row over the columns: each cell is hypergeometric
            // given the row's and the remaining columns' totals.
            let mut row_left = r;
            let mut pool = total_left;
            let (last, head) = col_left.split_last_mut().expect("two columns");
            for c in head {
                let x = hypergeometric(pool, *c, row_left, fact, rng);
                stat += fact.get(x);
                pool -= *c;
                *c -= x;
                row_left -= x;
            }
            stat += fact.get(row_left);
            *last -= row_left;
            total_left -= r;
        }
        // The last row takes what is left of the columns.
        stat += col_left.iter().map(|&c| fact.get(c)).sum::<f64>();
        stat >= self.threshold
    }
}

/// Number of successes in `draws` draws without replacement from `pool`
/// items of which `good` are successes.
fn hypergeometric(
    pool: i64,
    good: i64,
    draws: i64,
    fact: &LnFactorial,
    rng: &mut Xoshiro256pp,
) -> i64 {
    let bad = pool - good;
    let lo = (draws - bad).max(0);
    let hi = draws.min(good);
    if lo == hi {
        return lo;
    }
    let mode = ((i128::from(draws) + 1) * (i128::from(good) + 1) / (i128::from(pool) + 2)) as i64;
    let mode = mode.clamp(lo, hi);
    let p_mode = fact.ln_hypergeometric(mode, good, bad, draws).exp();

    // Inversion, visiting values outward from the mode.
    let mut u = rng.uniform() - p_mode;
    if u < 0. {
        return mode;
    }
    let (mut up, mut p_up) = (mode, p_mode);
    let (mut down, mut p_down) = (mode, p_mode);
    let f = |x: i64| x as f64;
    loop {
        let (can_up, can_down) = (up < hi, down > lo);
        if !(can_up || can_down) {
            // Rounding left `u` slightly above the total probability.
            return mode;
        }
        if can_up {
            p_up *= f(good - up) * f(draws - up) / (f(up + 1) * f(bad - draws + up + 1));
            up += 1;
            u -= p_up;
            if u < 0. {
                return up;
            }
        }
        if can_down {
            p_down *= f(down) * f(bad - draws + down) / (f(good - down + 1) * f(draws - down + 1));
            down -= 1;
            u -= p_down;
            if u < 0. {
                return down;
            }
        }
    }
}

/// `ln(n!)`: from a table up to [`LnFactorial::TABLE_MAX`], and from
/// Stirling's series above, so memory stays bounded for huge totals.
struct LnFactorial {
    table: LogFactorials,
    len: i64,
}

impl LnFactorial {
    /// Largest `n` kept in the table (8 MiB).
    const TABLE_MAX: i64 = 1 << 20;

    fn new(n: i64) -> Self {
        let len = n.min(Self::TABLE_MAX);
        LnFactorial {
            table: LogFactorials::new(len as i32),
            len,
        }
    }

    #[inline]
    fn get(&self, n: i64) -> f64 {
        if n <= self.len {
            self.table[n as i32]
        } else {
            let x = n as f64;
            (x + 0.5) * x.ln() - x + LN_SQRT_2PI + self.stirling_error(n)
        }
    }

    /// `ln P(X = x)` for X hypergeometric: `x` successes in `draws` draws
    /// from `good` successes and `bad` failures.
    fn ln_hypergeometric(&self, x: i64, good: i64, bad: i64, draws: i64) -> f64 {
        if good + bad <= self.len {
            // Small enough for the table, whose terms are small enough not
            // to cancel badly.
            let t = |n: i64| self.table[n as i32];
            t(good) + t(bad) + t(draws) + t(good + bad - draws)
                - t(good + bad)
                - t(x)
                - t(good - x)
                - t(draws - x)
                - t(bad - draws + x)
        } else {
            self.ln_hypergeometric_large(x, good, bad, draws)
        }
    }

    /// [`ln_hypergeometric`](Self::ln_hypergeometric) where `ln(n!)` is too
    /// large to subtract accurately: Loader's binomial-deviance form, as in
    /// R's `dhyper`.
    fn ln_hypergeometric_large(&self, x: i64, good: i64, bad: i64, draws: i64) -> f64 {
        let n = (good + bad) as f64;
        let p = draws as f64 / n;
        let q = (good + bad - draws) as f64 / n;
        self.ln_binomial(x, good, p, q) + self.ln_binomial(draws - x, bad, p, q)
            - self.ln_binomial(draws, good + bad, p, q)
    }

    /// `ln P(X = x)` for X binomial(n, p), `q = 1 - p` (Loader 2000, R's
    /// `dbinom_raw`).
    fn ln_binomial(&self, x: i64, n: i64, p: f64, q: f64) -> f64 {
        let (xf, nf) = (x as f64, n as f64);
        if x == 0 {
            return if n == 0 {
                0.
            } else if p < 0.1 {
                -bd0(nf, nf * q) - nf * p
            } else {
                nf * q.ln()
            };
        }
        if x == n {
            return if q < 0.1 {
                -bd0(nf, nf * p) - nf * q
            } else {
                nf * p.ln()
            };
        }
        let lc = self.stirling_error(n)
            - self.stirling_error(x)
            - self.stirling_error(n - x)
            - bd0(xf, nf * p)
            - bd0(nf - xf, nf * q);
        lc - 0.5 * (LN_2PI + xf.ln() + (-xf / nf).ln_1p())
    }

    /// `ln(n!) - ((n + 1/2) ln(n) - n + ln(sqrt(2 pi)))` for `n >= 1`.
    fn stirling_error(&self, n: i64) -> f64 {
        if n <= 15 {
            let x = n as f64;
            return self.table[n as i32] - ((x + 0.5) * x.ln() - x + LN_SQRT_2PI);
        }
        const S0: f64 = 1. / 12.;
        const S1: f64 = 1. / 360.;
        const S2: f64 = 1. / 1260.;
        const S3: f64 = 1. / 1680.;
        const S4: f64 = 1. / 1188.;
        let x = n as f64;
        let xx = x * x;
        if n > 500 {
            (S0 - S1 / xx) / x
        } else if n > 80 {
            (S0 - (S1 - S2 / xx) / xx) / x
        } else if n > 35 {
            (S0 - (S1 - (S2 - S3 / xx) / xx) / xx) / x
        } else {
            (S0 - (S1 - (S2 - (S3 - S4 / xx) / xx) / xx) / xx) / x
        }
    }
}

const LN_2PI: f64 = 1.837_877_066_409_345_5;
const LN_SQRT_2PI: f64 = 0.918_938_533_204_672_8;

/// The deviance term `x ln(x / np) + np - x`, accurately when `x` is close
/// to `np` (Loader 2000, R's `bd0`).
fn bd0(x: f64, np: f64) -> f64 {
    if (x - np).abs() < 0.1 * (x + np) {
        let v = (x - np) / (x + np);
        let mut s = (x - np) * v;
        let mut ej = 2. * x * v;
        let v2 = v * v;
        for j in 1..1000 {
            ej *= v2;
            let s1 = s + ej / f64::from(2 * j + 1);
            if s1 == s {
                return s1;
            }
            s = s1;
        }
        s
    } else {
        x * (x / np).ln() + np - x
    }
}

/// xoshiro256++ (Blackman & Vigna).
struct Xoshiro256pp([u64; 4]);

impl Xoshiro256pp {
    /// The generator of stream `stream` for `seed`.
    fn new(seed: u64, stream: u64) -> Self {
        let mut state = splitmix64(seed ^ splitmix64(stream));
        let mut next = || {
            state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            splitmix64(state)
        };
        Xoshiro256pp([next(), next(), next(), next()])
    }

    fn next_u64(&mut self) -> u64 {
        let s = &mut self.0;
        let result = s[0].wrapping_add(s[3]).rotate_left(23).wrapping_add(s[0]);
        let t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = s[3].rotate_left(45);
        result
    }

    /// Uniform in [0, 1).
    fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1. / (1u64 << 53) as f64)
    }
}

/// The SplitMix64 output function.
fn splitmix64(x: u64) -> u64 {
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{fisher_exact, Options};

    fn simulate(rows: &[&[i64]], samples: u64, seed: u64, threads: usize) -> Simulated {
        let cells: Vec<i64> = rows.concat();
        let table = Table::from_row_major(&cells, rows.len(), rows[0].len());
        let simulation = Simulation {
            samples,
            seed: Some(seed),
            threads: NonZeroUsize::new(threads),
        };
        fisher_simulated(&table, &simulation).unwrap()
    }

    fn exact(rows: &[&[i64]]) -> f64 {
        let cells: Vec<i64> = rows.concat();
        let table = Table::from_row_major(&cells, rows.len(), rows[0].len());
        fisher_exact(&table, &Options::default()).unwrap().p_value
    }

    #[test]
    fn hypergeometric_distribution() {
        // pool 20, 7 good, 9 draws: compare frequencies with the pmf, with
        // the log-factorial table and with the large-number formulas.
        check_hypergeometric(&LnFactorial::new(20));
        check_hypergeometric(&small_table(15));
    }

    fn check_hypergeometric(fact: &LnFactorial) {
        let mut rng = Xoshiro256pp::new(1, 0);
        let n = 400_000;
        let mut counts = [0u32; 8];
        for _ in 0..n {
            counts[hypergeometric(20, 7, 9, fact, &mut rng) as usize] += 1;
        }
        let exact = LogFactorials::new(20);
        let choose = |n: i32, k: i32| (exact[n] - exact[k] - exact[n - k]).exp();
        for (x, &c) in counts.iter().enumerate() {
            let x = x as i32;
            let p = choose(7, x) * choose(13, 9 - x) / choose(20, 9);
            let sd = (p * (1. - p) / n as f64).sqrt();
            let freq = f64::from(c) / n as f64;
            assert!((freq - p).abs() <= 5. * sd + 1e-12, "x={x}: {freq} vs {p}");
        }
    }

    /// A table of `len` entries, using the large-number formulas above.
    fn small_table(len: i64) -> LnFactorial {
        LnFactorial {
            table: LogFactorials::new(len as i32),
            len,
        }
    }

    #[test]
    fn stirling_matches_table() {
        let full = LogFactorials::new(20_000);
        let fact = small_table(100);
        for n in (101..20_000).step_by(7) {
            let (a, b) = (fact.get(n), full[n as i32]);
            // (The table accumulates rounding errors, ~1e-13 at n = 10^4.)
            assert!((a - b).abs() <= 1e-12 * b, "{n}: {a} vs {b}");
        }
        // ln(353!) to 20 digits (mpmath).
        assert!((fact.get(353) - 1_721.715_632_730_827_7).abs() <= 1e-15 * 1721.7);
    }

    #[test]
    fn loader_matches_table() {
        let fact = LnFactorial::new(3000);
        for (good, bad, draws) in [
            (7, 13, 9),
            (40, 2900, 60),
            (1500, 1500, 1500),
            (3, 2000, 1000),
        ] {
            let lo = (draws - bad).max(0);
            for x in lo..=draws.min(good) {
                let a = fact.ln_hypergeometric(x, good, bad, draws);
                let b = fact.ln_hypergeometric_large(x, good, bad, draws);
                // (The table's terms cancel: it is the less accurate one.)
                assert!((a - b).abs() <= 1e-9 * a.abs().max(1.), "{x}: {a} vs {b}");
            }
        }
        // ln P(X = 730), X ~ hypergeometric(1500, 1500, 1500), to 20 digits.
        let v = fact.ln_hypergeometric_large(730, 1500, 1500, 1500);
        assert!((v - -4.602_160_030_483_242).abs() <= 1e-14, "{v}");
    }

    #[test]
    fn huge_counts() {
        // 2 x 2 with total 4e9: x11 ~ N(m, m / 4) approximately, so a cell
        // two standard deviations off the mean gives p ~ 0.0455.
        let m: i64 = 1_000_000_000;
        let d = 31_623; // 2 * sqrt(m) / 2
        let p = simulate(&[&[m + d, m - d], &[m - d, m + d]], 20_000, 5, 0).p_value;
        assert!((p - 0.0455).abs() < 0.006, "{p}");
    }

    #[test]
    fn agrees_with_exact() {
        let tables: &[&[&[i64]]] = &[
            &[&[5, 1, 5], &[1, 1, 5]],
            &[&[10, 6, 7], &[5, 7, 9]],
            &[&[2, 1, 3, 0], &[1, 2, 1, 4], &[0, 3, 1, 2]],
            &[&[8, 2], &[1, 5]],
            &[&[20, 3, 0, 4], &[2, 15, 6, 1], &[1, 4, 12, 3]],
        ];
        for t in tables {
            let p = exact(t);
            let sim = simulate(t, 200_000, 7, 0);
            let tol = 5. * (p * (1. - p) / 200_000.).sqrt() + 1e-5;
            assert!(
                (sim.p_value - p).abs() <= tol,
                "{t:?}: {} vs {p}",
                sim.p_value
            );
        }
    }

    #[test]
    fn independent_of_threads() {
        let t: &[&[i64]] = &[&[10, 6, 7], &[5, 7, 9], &[3, 8, 2]];
        let one = simulate(t, 10_000, 42, 1);
        assert_eq!(one, simulate(t, 10_000, 42, 3));
        assert_eq!(one, simulate(t, 10_000, 42, 16));
        assert_ne!(one, simulate(t, 10_000, 43, 1));
    }

    #[test]
    fn edge_cases() {
        let one_row = simulate(&[&[3, 0, 4]], 100, 1, 0);
        assert_eq!(one_row.p_value, 1.);
        // Zero rows and columns are ignored.
        let a = simulate(&[&[8, 0, 2], &[0, 0, 0], &[1, 0, 5]], 5000, 3, 0);
        let b = simulate(&[&[8, 2], &[1, 5]], 5000, 3, 0);
        assert_eq!(a, b);
        assert_eq!(simulate(&[&[8, 2], &[1, 5]], 0, 3, 0).p_value, 1.);
        let cells = [1, -1];
        let table = Table::from_row_major(&cells, 1, 2);
        assert_eq!(
            fisher_simulated(&table, &Simulation::default()),
            Err(Error::NegativeEntry)
        );
    }
}
