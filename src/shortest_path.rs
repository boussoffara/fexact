//! Shortest path through the remaining network (C: `f4xact`, `f8xact`,
//! `f11act`).
//!
//! Finds (an upper bound on) the least probable completion of a subtable,
//! i.e. the maximum of `sum ln(cell!)`, by depth-first search over greedy
//! fillings: each step puts `min(row, column)` into one cell, which zeroes a
//! row or column (or both) and shrinks the table. The search stops early once
//! the bound cannot beat `dspt`.

use crate::math::LogFactorials;
use crate::TOL;

/// Search state of one level (C: `nrstk`, `ncstk`, `ystk`, `lstk`, `mstk`,
/// `nstk`).
#[derive(Debug, Clone, Copy, Default)]
struct Frame {
    /// Rows and columns left at this level.
    nrows: usize,
    ncols: usize,
    /// `sum ln(cell!)` of the cells filled to reach it.
    y: f64,
    /// The candidate being tried, and the end of the candidates.
    cand: usize,
    end: usize,
    /// Whether candidates are rows (paired with the first column) or
    /// columns (paired with the first row).
    vary_rows: bool,
}

/// Reusable work arrays for [`ShortestPath::compute`].
pub(crate) struct ShortestPath {
    /// Row and column totals of each level, in descending order (C:
    /// `irstk`, `icstk`).
    row_stack: Vec<i32>,
    col_stack: Vec<i32>,
    frames: Vec<Frame>,
}

impl ShortestPath {
    /// Work arrays for subtables with at most `levels` (rows + columns + 1)
    /// levels and `dim` rows or columns.
    pub(crate) fn new(levels: usize, dim: usize) -> Self {
        ShortestPath {
            row_stack: vec![0; levels * dim],
            col_stack: vec![0; levels * dim],
            frames: vec![Frame::default(); levels],
        }
    }

    /// Shortest path length for a table with sorted row totals `irow` and
    /// sorted column totals `icol`, as an offset from `dspt`; returns
    /// `-dspt` when no path can be shorter than `dspt`.
    pub(crate) fn compute(
        &mut self,
        irow: &[i32],
        icol: &[i32],
        dspt: f64,
        fact: &LogFactorials,
    ) -> f64 {
        let (nrow, ncol) = (irow.len(), icol.len());
        let neg_sum = |xs: &[i32]| xs.iter().fold(0., |sp, &x| sp - fact[x]);
        match (nrow, ncol) {
            (1, _) => return neg_sum(icol),
            (_, 1) => return neg_sum(irow),
            (2, 2) => {
                // Put the smaller of the larger totals in one cell.
                let (a, b) = (irow[1], icol[1]);
                return if a <= b {
                    -(fact[a] + fact[b] + fact[b - a])
                } else {
                    -(fact[b] + fact[a] + fact[a - b])
                };
            }
            _ => {}
        }

        let ShortestPath {
            row_stack,
            col_stack,
            frames,
        } = self;
        // Totals of level `s`: rows from s * nrow, columns from s * ncol.
        for (dst, &x) in row_stack.iter_mut().zip(irow.iter().rev()) {
            *dst = x;
        }
        for (dst, &x) in col_stack.iter_mut().zip(icol.iter().rev()) {
            *dst = x;
        }

        let mut s = 0; // current level
        let (mut nro, mut nco) = (nrow, ncol);
        let mut y = 0.;
        frames[0] = Frame {
            nrows: nro,
            ncols: nco,
            y,
            ..Frame::default()
        };
        let mut cand = 0;
        let mut amx = 0.; // the largest complete y so far
        let mut sp = dspt;

        // `end` is deliberately not reloaded when resuming a level (as in
        // the Fortran original), so it lives outside the loops.
        let (mut end, mut vary_rows) = (0, true);
        let mut resume = false;

        loop {
            // Descend greedily until one row or column is left.
            loop {
                if !resume {
                    // Vary whichever side's largest total is smaller (rows on
                    // a tie).
                    let r1 = row_stack[s * nrow];
                    let c1 = col_stack[s * ncol];
                    (end, vary_rows) = if r1 > c1 {
                        if nro >= nco {
                            (nco - 1, false)
                        } else {
                            (nro, true)
                        }
                    } else if r1 < c1 {
                        if nro <= nco {
                            (nro - 1, true)
                        } else {
                            (nco, false)
                        }
                    } else if nro <= nco {
                        (nro - 1, true)
                    } else {
                        (nco - 1, false)
                    };
                }
                resume = false;

                // Fill cell (i, j) with min(row i, column j).
                let (i, j) = if vary_rows { (cand, 0) } else { (0, cand) };
                let (rows, next_rows) = row_stack.split_at_mut((s + 1) * nrow);
                let (cols, next_cols) = col_stack.split_at_mut((s + 1) * ncol);
                let (rows, cols) = (&rows[s * nrow..], &cols[s * ncol..]);
                let (irt, ict) = (rows[i], cols[j]);
                y += fact[irt.min(ict)];
                if irt == ict {
                    nro -= 1;
                    nco -= 1;
                    remove_at(rows, i, nro, next_rows);
                    remove_at(cols, j, nco, next_cols);
                } else if irt > ict {
                    nco -= 1;
                    remove_at(cols, j, nco, next_cols);
                    replace_at(rows, i, irt - ict, nro, next_rows);
                } else {
                    nro -= 1;
                    remove_at(rows, i, nro, next_rows);
                    replace_at(cols, j, ict - irt, nco, next_cols);
                }

                if nro == 1 {
                    y = next_cols[..nco].iter().fold(y, |y, &c| y + fact[c]);
                    break;
                }
                if nco == 1 {
                    y = next_rows[..nro].iter().fold(y, |y, &r| y + fact[r]);
                    break;
                }

                frames[s].cand = cand;
                frames[s].end = end;
                frames[s].vary_rows = vary_rows;
                s += 1;
                frames[s].nrows = nro;
                frames[s].ncols = nco;
                frames[s].y = y;
                cand = 0;
            }

            if y > amx {
                amx = y;
                if sp - amx <= TOL {
                    return -dspt;
                }
            }

            // Back up to the deepest earlier level with an untried candidate
            // that differs from the one before it.
            'backtrack: loop {
                let Some(prev) = s.checked_sub(1) else {
                    sp -= amx;
                    return if sp - amx <= TOL { -dspt } else { sp - dspt };
                };
                s = prev;
                let f = frames[s];
                cand = f.cand + 1;
                while cand < f.end {
                    (vary_rows, nro, nco, y) = (f.vary_rows, f.nrows, f.ncols, f.y);
                    let distinct = if vary_rows {
                        row_stack[s * nrow + cand] < row_stack[s * nrow + cand - 1]
                    } else {
                        col_stack[s * ncol + cand] < col_stack[s * ncol + cand - 1]
                    };
                    if distinct {
                        resume = true;
                        break 'backtrack;
                    }
                    cand += 1;
                }
            }
        }
    }
}

/// Writes `src[..=len]` without the element at `at` to `dst[..len]`
/// (C: `f11act`).
///
/// (Element-wise loops: the vectors are short, and `copy_from_slice` would
/// call `memcpy`.)
fn remove_at(src: &[i32], at: usize, len: usize, dst: &mut [i32]) {
    for (d, &x) in dst[..len]
        .iter_mut()
        .zip(src[..at].iter().chain(&src[at + 1..=len]))
    {
        *d = x;
    }
}

/// Writes the descending `src[..len]` with the element at `at` replaced by
/// the smaller `value` to `dst[..len]`, keeping it descending (C: `f8xact`).
fn replace_at(src: &[i32], at: usize, value: i32, len: usize, dst: &mut [i32]) {
    let (src, dst) = (&src[..len], &mut dst[..len]);
    // Shift the larger totals after `at` up one place, then insert.
    let mut i = at;
    while i + 1 < len && value < src[i + 1] {
        dst[i] = src[i + 1];
        i += 1;
    }
    dst[i] = value;
    for k in (0..at).chain(i + 1..len) {
        dst[k] = src[k];
    }
}
