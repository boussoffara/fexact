//! Fisher's exact test for r x c contingency tables.
//!
//! ```
//! use fexact::{fisher_exact, Options, Table};
//!
//! let counts = [5, 1, 5,
//!               1, 1, 5];
//! let table = Table::from_row_major(&counts, 2, 3);
//! let out = fisher_exact(&table, &Options::default()).unwrap();
//! assert!((out.p_value - 0.4004524886877815).abs() < 1e-12);
//! ```
//!
//! A port of `legacy/fexact.c`, R's version of ACM TOMS Algorithm 643
//! (Mehta & Patel's network algorithm, 1993), producing bit-identical
//! p-values and the same errors (see `tests/golden.rs`), except where the C
//! code's 32-bit integers overflow: node keys and path counts are 64-bit
//! here, so tables the C code rejects (error 501) or silently gets wrong
//! (wide tables, whose path counts wrap around) are computed correctly.
//! Modules map to the C routines:
//!
//! | module            | C                                     | role                        |
//! |-------------------|---------------------------------------|-----------------------------|
//! | [`compat`]        | `fexact`                              | the C interface             |
//! | `workspace`       | `fexact`, `iwork`                     | workspace budget            |
//! | `network`         | `f2xact`                              | the network algorithm       |
//! | `past_paths`      | `f5xact`, `f6xact`                    | nodes and past path lengths |
//! | `nodes`           | `f7xact`                              | enumerating daughters       |
//! | `longest_path`    | `f3xact`, `f10act`                    | longest-path bound          |
//! | `shortest_path`   | `f4xact`, `f8xact`, `f11act`          | shortest-path bound         |
//! | `math`, `sort`    | `f9xact`, `gammds`, `alogam`, `isort` | numeric helpers             |
//! | `simulate`        | (R's `simulate.p.value`)              | Monte Carlo p-values        |
//!
//! Where the C code's behaviour is otherwise undefined (reading out of
//! bounds, dividing by zero), this port returns an error or, for one benign
//! out-of-bounds read, its practical value; see `math::LogFactorials` and
//! `workspace`.
//!
//! The optional `python` feature adds the `fexact` Python extension module.

mod error;
mod longest_path;
mod math;
mod network;
mod nodes;
mod past_paths;
#[cfg(feature = "python")]
mod python;
mod shortest_path;
mod simulate;
mod sort;
mod workspace;

use std::num::NonZeroUsize;
use std::sync::atomic::{AtomicUsize, Ordering};

pub use crate::error::Error;
pub use crate::simulate::{fisher_simulated, Simulated, Simulation};

use crate::workspace::Layout;

/// Tolerance for treating path lengths as equal: the square root of the
/// smallest relative spacing (single precision, as in the original).
pub(crate) const TOL: f64 = 3.45254e-7;

/// Encoded vector of margins, the key of the hash tables of nodes. (64-bit:
/// the C code's `int` keys overflow on wide tables.)
pub(crate) type Key = i64;

/// Hash slot of `key` in a table of `n` slots.
pub(crate) fn slot_of(key: Key, n: usize) -> usize {
    (key as u64 % n as u64) as usize
}

/// A read-only view of an `nrows` x `ncols` table of counts.
#[derive(Debug, Clone, Copy)]
pub struct Table<'a> {
    cells: &'a [i64],
    nrows: usize,
    ncols: usize,
    row_stride: usize,
    col_stride: usize,
}

impl<'a> Table<'a> {
    /// A table whose rows are stored one after another.
    ///
    /// # Panics
    ///
    /// If `cells.len() != nrows * ncols`.
    pub fn from_row_major(cells: &'a [i64], nrows: usize, ncols: usize) -> Self {
        assert_eq!(cells.len(), nrows * ncols, "cells do not match the shape");
        Table {
            cells,
            nrows,
            ncols,
            row_stride: ncols,
            col_stride: 1,
        }
    }

    /// A table whose columns are stored one after another.
    ///
    /// # Panics
    ///
    /// If `cells.len() != nrows * ncols`.
    pub fn from_column_major(cells: &'a [i64], nrows: usize, ncols: usize) -> Self {
        assert_eq!(cells.len(), nrows * ncols, "cells do not match the shape");
        Table {
            cells,
            nrows,
            ncols,
            row_stride: 1,
            col_stride: nrows,
        }
    }

    pub fn nrows(&self) -> usize {
        self.nrows
    }

    pub fn ncols(&self) -> usize {
        self.ncols
    }

    /// The count in row `i`, column `j`.
    ///
    /// # Panics
    ///
    /// If the cell is out of bounds.
    pub fn get(&self, i: usize, j: usize) -> i64 {
        assert!(i < self.nrows && j < self.ncols, "cell out of bounds");
        self.cells[i * self.row_stride + j * self.col_stride]
    }

    /// [`get`](Self::get) as an `i32`, for the exact test once it has
    /// checked that the table total fits.
    pub(crate) fn cell(&self, i: usize, j: usize) -> i32 {
        self.get(i, j) as i32
    }
}

/// Parameters of the hybrid approximation: a subtable's completions are
/// approximated by a chi-squared tail when at least `percnt` percent of its
/// cells have expected count `expect` or more, and none has less than
/// `emin`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hybrid {
    pub expect: f64,
    pub percnt: f64,
    pub emin: f64,
}

impl Hybrid {
    /// Cochran's conditions: expected counts of at least 5 in 80% of the
    /// cells and at least 1 everywhere.
    pub const COCHRAN: Hybrid = Hybrid {
        expect: 5.,
        percnt: 80.,
        emin: 1.,
    };
}

/// Settings for [`fisher_exact`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Options {
    /// Workspace size in 4-byte units, which bounds the problem size the
    /// algorithm can handle (and its memory use). Default 200 000.
    pub workspace: u32,
    /// Use the hybrid approximation instead of the exact test.
    /// Default `None`.
    pub hybrid: Option<Hybrid>,
    /// Ratio of past-path entries to node slots in the workspace split.
    /// Default 30.
    pub mult: u32,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            workspace: 200_000,
            hybrid: None,
            mult: 30,
        }
    }
}

/// Result of [`fisher_exact`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Output {
    /// Probability of a table at most as likely as the observed one, given
    /// the margins.
    pub p_value: f64,
    /// Probability of the observed table given the margins.
    pub table_probability: f64,
}

/// Fisher's exact test (or its hybrid approximation) for `table`.
pub fn fisher_exact(table: &Table, options: &Options) -> Result<Output, Error> {
    let mut ntot: i64 = 0;
    for j in 0..table.ncols {
        for i in 0..table.nrows {
            let x = table.get(i, j);
            if x < 0 {
                return Err(Error::NegativeEntry);
            }
            ntot = ntot.checked_add(x).ok_or(Error::OutOfWorkspace)?;
        }
    }
    if ntot == 0 {
        return Err(Error::AllZero);
    }
    // The log-factorial table alone needs 2 * (ntot + 1) ints of workspace,
    // so every cell fits an i32 from here on.
    let ntot = i32::try_from(ntot)
        .ok()
        .filter(|&n| n < i32::MAX)
        .ok_or(Error::OutOfWorkspace)?;

    let layout = Layout::new(
        options.workspace.into(),
        options.mult.into(),
        table.nrows,
        table.ncols,
        ntot,
    )?;
    network::network(table, options.hybrid.as_ref(), &layout)
}

/// [`fisher_exact`] for each of `tables`, on `threads` threads (default:
/// all available cores). Results are in the order of `tables`.
pub fn fisher_exact_many(
    tables: &[Table],
    options: &Options,
    threads: Option<NonZeroUsize>,
) -> Vec<Result<Output, Error>> {
    let threads = threads
        .or_else(|| std::thread::available_parallelism().ok())
        .map_or(1, NonZeroUsize::get)
        .min(tables.len());
    let next = AtomicUsize::new(0);
    let mut results = vec![None; tables.len()];
    std::thread::scope(|s| {
        let workers: Vec<_> = (0..threads)
            .map(|_| {
                s.spawn(|| {
                    // Tables vary wildly in cost: hand them out one by one.
                    let mut done = Vec::new();
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        let Some(table) = tables.get(i) else {
                            return done;
                        };
                        done.push((i, fisher_exact(table, options)));
                    }
                })
            })
            .collect();
        for worker in workers {
            let done = worker
                .join()
                .unwrap_or_else(|panic| std::panic::resume_unwind(panic));
            for (i, result) in done {
                results[i] = Some(result);
            }
        }
    });
    results
        .into_iter()
        .map(|r| r.expect("every table is processed"))
        .collect()
}

/// The C library's interface.
pub mod compat {
    use super::*;

    /// Computes Fisher's exact test (or the hybrid approximation when
    /// `expect > 0`) for an `nrow` x `ncol` table stored column-major with
    /// leading dimension `ldtabl`, i.e. cell (i, j) is
    /// `table[i + j * ldtabl]`; `workspace` and `mult` as in R (200000, 30).
    ///
    /// # Panics
    ///
    /// If `nrow`, `ncol` or `ldtabl` is negative, or `table` holds fewer
    /// than `ldtabl * ncol` cells.
    #[allow(clippy::too_many_arguments)]
    pub fn fexact(
        nrow: i32,
        ncol: i32,
        table: &[i32],
        ldtabl: i32,
        expect: f64,
        percnt: f64,
        emin: f64,
        workspace: i32,
        mult: i32,
    ) -> Result<Output, Error> {
        // C: malloc() of a negative size.
        if 2 * (workspace / 2) < 0 {
            return Err(Error::AllocationFailed);
        }
        if nrow > ldtabl {
            return Err(Error::LeadingDimension);
        }
        let dim = |x: i32| usize::try_from(x).expect("dimensions must not be negative");
        let (nrows, ncols, ld) = (dim(nrow), dim(ncol), dim(ldtabl));
        let cells: Vec<i64> = table[..ld * ncols].iter().map(|&x| x.into()).collect();
        let table = Table {
            cells: &cells,
            nrows,
            ncols,
            row_stride: 1,
            col_stride: ld,
        };
        let options = Options {
            // -1 behaves like 0 in C (the workspace is rounded down to even).
            workspace: workspace.max(0) as u32,
            hybrid: (expect > 0.).then_some(Hybrid {
                expect,
                percnt,
                emin,
            }),
            mult: u32::try_from(mult).expect("mult must not be negative"),
        };
        fisher_exact(&table, &options)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(rows: &[&[i64]], workspace: u32, hybrid: bool) -> Result<f64, Error> {
        let cells: Vec<i64> = rows.concat();
        let table = Table::from_row_major(&cells, rows.len(), rows[0].len());
        let options = Options {
            workspace,
            hybrid: hybrid.then_some(Hybrid::COCHRAN),
            ..Options::default()
        };
        fisher_exact(&table, &options).map(|r| r.p_value)
    }

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-12 * b.abs().max(1e-300)
    }

    #[test]
    fn small_2x3() {
        let p = run(&[&[5, 1, 5], &[1, 1, 5]], 200_000, false).unwrap();
        assert!(close(p, 0.4004524886877815), "{p}");
    }

    #[test]
    fn transpose_invariant() {
        let a = run(
            &[&[2, 1, 3, 0], &[1, 2, 1, 4], &[0, 3, 1, 2]],
            200_000,
            false,
        )
        .unwrap();
        let b = run(
            &[&[2, 1, 0], &[1, 2, 3], &[3, 1, 1], &[0, 4, 2]],
            200_000,
            false,
        )
        .unwrap();
        assert!(close(a, b), "{a} vs {b}");
    }

    #[test]
    fn row_and_column_major_agree() {
        let cells = [3, 1, 4, 1, 5, 9];
        let a = fisher_exact(&Table::from_row_major(&cells, 2, 3), &Options::default());
        let transposed = [3, 1, 1, 5, 4, 9];
        let b = fisher_exact(
            &Table::from_column_major(&transposed, 2, 3),
            &Options::default(),
        );
        assert_eq!(a, b);
    }

    #[test]
    fn exact_vs_hybrid() {
        let t: &[&[i64]] = &[&[10, 6, 7], &[5, 7, 9]];
        let exact = run(t, 10_000_000, false).unwrap();
        let hybrid = run(t, 10_000_000, true).unwrap();
        assert!((exact - 0.4158).abs() < 1e-3, "{exact}");
        assert!((hybrid - 0.3962).abs() < 1e-3, "{hybrid}");
    }

    #[test]
    fn errors() {
        assert_eq!(
            run(&[&[-1, 6], &[5, 7]], 200_000, false),
            Err(Error::NegativeEntry)
        );
        assert_eq!(
            run(&[&[0, 0], &[0, 0]], 200_000, false),
            Err(Error::AllZero)
        );
        assert_eq!(
            run(&[&[1], &[2]], 200_000, false),
            Err(Error::TooFewColumns)
        );
        assert_eq!(
            run(&[&[10, 6, 7], &[5, 7, 9]], 8, false),
            Err(Error::OutOfWorkspace)
        );
    }

    #[test]
    fn compat_errors() {
        let t = [10, 5, 6, 7, 7, 9];
        let run = |ws| compat::fexact(2, 3, &t, 2, -1., 100., 0., ws, 30).map(|r| r.p_value);
        assert_eq!(run(-5), Err(Error::AllocationFailed));
        assert_eq!(
            compat::fexact(3, 2, &t, 2, -1., 100., 0., 200_000, 30),
            Err(Error::LeadingDimension)
        );
        assert!((run(200_000).unwrap() - 0.4158).abs() < 1e-3);
        // C treats -1 like 0: an empty workspace.
        assert_eq!(run(-1), Err(Error::OutOfWorkspace));
    }

    #[test]
    fn many_matches_single() {
        let tables: Vec<Vec<i64>> = (1..40)
            .map(|i| vec![i % 7, 3, i % 5 + 1, 2, 4, i % 3])
            .collect();
        let views: Vec<Table> = tables
            .iter()
            .map(|c| Table::from_row_major(c, 2, 3))
            .collect();
        let options = Options::default();
        let many = fisher_exact_many(&views, &options, NonZeroUsize::new(4));
        let single: Vec<_> = views.iter().map(|t| fisher_exact(t, &options)).collect();
        assert_eq!(many, single);
    }
}
