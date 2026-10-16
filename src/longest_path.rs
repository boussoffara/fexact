//! Longest path through the remaining network (C: `f3xact`, `f10act`).
//!
//! For a subtable with given row and column totals, finds the completion
//! minimising `sum(ln cell!)`, i.e. the most probable table; the negated
//! minimum bounds the length of every remaining path, so past paths that
//! cannot exceed the observed length can be counted wholesale.
//!
//! The search fills one row at a time. Within a row it enumerates cell
//! values column by column between bounds derived from the expected counts;
//! each resulting vector of remaining column totals is a node of the next
//! row, kept with its best value so far in a [`MarginStack`].

use crate::error::Error;
use crate::math::{c_fmin, LogFactorials};
use crate::sort::sort_ascending;
use crate::{slot_of, Key, TOL};

/// Hash set of column-total vectors (encoded as keys) with their best
/// partial values, in two halves: the row being expanded (popped from) and
/// the next row (pushed to). C: `ist` (keys), `stv` (values), `itc` (order).
struct MarginStack {
    keys: Vec<Key>,
    values: Vec<f64>,
    /// Slots of the entries of each half, in insertion order.
    order: Vec<usize>,
    /// Slots per half.
    half: usize,
    /// Entries pushed to the next row (C: `nst`).
    pushed: usize,
    /// Entries of the current row not yet popped (C: `nitc`).
    pending: usize,
    /// Offsets of the halves being pushed to and popped from.
    push_half: usize,
    pop_half: usize,
}

const EMPTY: Key = -1;

impl MarginStack {
    fn new(half: usize) -> Self {
        MarginStack {
            keys: vec![EMPTY; 2 * half],
            values: vec![0.; 2 * half],
            order: vec![0; 2 * half],
            half,
            pushed: 0,
            pending: 0,
            push_half: 0,
            pop_half: half,
        }
    }

    fn reset(&mut self) {
        // A completed search pops every entry it pushed, so the slots are
        // already empty (after an error the computation is abandoned).
        debug_assert!(self.keys.iter().all(|&k| k == EMPTY));
        self.pushed = 0;
        self.pending = 0;
        self.push_half = 0;
        self.pop_half = self.half;
    }

    /// Adds `key` to the next row with value `v`, or lowers its value.
    fn push(&mut self, key: Key, v: f64) -> Result<(), Error> {
        let start = slot_of(key, self.half);
        for slot in (start..self.half).chain(0..start) {
            let i = self.push_half + slot;
            if self.keys[i] == EMPTY {
                self.keys[i] = key;
                self.values[i] = v;
                self.order[self.push_half + self.pushed] = slot;
                self.pushed += 1;
                return Ok(());
            } else if self.keys[i] == key {
                self.values[i] = c_fmin(v, self.values[i]);
                return Ok(());
            }
        }
        Err(Error::MarginStackFull)
    }

    /// Removes the most recently pushed entry of the current row.
    fn pop(&mut self) -> Option<(Key, f64)> {
        self.pending = self.pending.checked_sub(1)?;
        let i = self.pop_half + self.order[self.pop_half + self.pending];
        let key = std::mem::replace(&mut self.keys[i], EMPTY);
        Some((key, self.values[i]))
    }

    /// Makes the entries pushed so far the current row.
    fn next_level(&mut self) {
        self.pending = self.pushed;
        self.pushed = 0;
        self.pop_half = self.push_half;
        self.push_half = self.half - self.push_half;
    }
}

/// The search position within one column of the row being filled.
#[derive(Debug, Clone, Copy, Default)]
struct Cell {
    /// Value of the cell (C: `lb`).
    value: i32,
    /// Values left to try, counting the current one (C: `nu`).
    tries: i32,
    /// Row total left after this cell (C: `nr`).
    row_left: i32,
    /// Total of the columns after this one (C: `nt`).
    cols_after: i32,
}

/// Reusable work arrays for [`LongestPath::compute`].
pub(crate) struct LongestPath {
    stack: MarginStack,
    /// Row totals; rows `first..` are still to be filled (C: `iro`).
    rows: Vec<i32>,
    /// Column totals of the current node (C: `ico`).
    cols: Vec<i32>,
    cells: Vec<Cell>,
    /// `prefix[j]`: sum of `ln(value!)` over the row's first `j` cells
    /// (C: `alen`).
    prefix: Vec<f64>,
    /// Scratch for building keys, and for [`special_table`].
    key: Vec<i32>,
    scratch: SpecialScratch,
}

impl LongestPath {
    /// Work arrays for subtables with at most `dim` rows and columns, with a
    /// stack of `2 * half` entries.
    pub(crate) fn new(dim: usize, half: usize) -> Self {
        LongestPath {
            stack: MarginStack::new(half),
            rows: vec![0; dim],
            cols: vec![0; dim],
            cells: vec![Cell::default(); dim],
            prefix: vec![0.; dim + 1],
            key: Vec::with_capacity(dim),
            scratch: SpecialScratch::default(),
        }
    }

    /// Longest path length for a table with sorted row totals `irow` and
    /// sorted column totals `icol`, both summing to `ntot`.
    pub(crate) fn compute(
        &mut self,
        irow: &[i32],
        icol: &[i32],
        ntot: i32,
        fact: &LogFactorials,
    ) -> Result<f64, Error> {
        if let Some(v) = trivial(irow, icol, ntot, fact) {
            return Ok(v);
        }
        if let Some(v) = closed_form(irow, icol, 0., fact, &mut self.scratch) {
            return Ok(-v);
        }
        self.search(irow, icol, ntot, fact)
    }

    /// The row-by-row search (the general case of `f3xact`).
    fn search(
        &mut self,
        irow: &[i32],
        icol: &[i32],
        ntot: i32,
        fact: &LogFactorials,
    ) -> Result<f64, Error> {
        // Fill along the longer dimension: no more columns than rows.
        let (rows_in, cols_in) = if irow.len() >= icol.len() {
            (irow, icol)
        } else {
            (icol, irow)
        };
        let ncols = cols_in.len();
        let LongestPath {
            stack,
            rows,
            cols,
            cells,
            prefix,
            key,
            scratch,
        } = self;
        let (rows, cols, cells) = (
            &mut rows[..rows_in.len()],
            &mut cols[..ncols],
            &mut cells[..ncols],
        );
        rows.copy_from_slice(rows_in);
        cols.copy_from_slice(cols_in);
        prefix[..=icol.len()].fill(0.);
        stack.reset();

        let mut first = 0; // the row being filled
        let mut rows_left = rows.len();
        let mut total = ntot; // of rows first..
        set_cols_after(cells, cols, total);
        // Radix for encoding column totals (the largest is the last).
        let radix = Key::from(cols[ncols - 1]) + 1;
        let mut base = 0.; // value of the node being expanded (C: `val`)
        let mut best = 1e100; // min over completions of sum ln(cell!) (C: `vmn`)
        let mut j = 0; // the column being varied

        enum Step {
            /// Start filling row `first` for the node in `cols`.
            NewNode,
            /// Try the next value of cell `j`, backtracking as needed.
            NextValue,
            /// Take the next node of this row, or move to the next row.
            PopNode,
        }
        let mut step = Step::NewNode;

        loop {
            match step {
                Step::NewNode => {
                    // The first cell starts one below its bound: NextValue
                    // increments it before use.
                    j = 0;
                    let (nrt, nct) = (rows[first], cols[0]);
                    let others = rows_left as i32 - 1;
                    let c = &mut cells[0];
                    c.value = ((f64::from(nrt) + 1.) * f64::from(nct + 1)
                        / f64::from(total + others * (ncols as i32 - 1) + 1)
                        - TOL) as i32
                        - 1;
                    c.tries = ((f64::from(nrt) + f64::from(ncols as i32 - 1))
                        * f64::from(nct + others)
                        / f64::from(total + others + ncols as i32 - 1))
                        as i32
                        - c.value
                        + 1;
                    c.row_left = nrt - c.value;
                    step = Step::NextValue;
                }
                Step::NextValue => {
                    cells[j].tries -= 1;
                    if cells[j].tries == 0 {
                        if j == 0 {
                            step = Step::PopNode;
                        } else {
                            j -= 1;
                        }
                        continue;
                    }
                    cells[j].value += 1;
                    cells[j].row_left -= 1;

                    // Set the following cells to their lower bounds.
                    let others = rows_left as i32 - 1;
                    loop {
                        prefix[j + 1] = prefix[j] + fact[cells[j].value];
                        if j + 2 >= ncols {
                            break;
                        }
                        let Cell {
                            row_left: nrt,
                            cols_after: nn1,
                            ..
                        } = cells[j];
                        j += 1;
                        let nc1 = (ncols - j - 1) as i32;
                        let nct = cols[j];
                        let c = &mut cells[j];
                        c.value = ((f64::from(nrt) + 1.) * f64::from(nct + 1)
                            / f64::from(nn1 + others * nc1 + 1)
                            - TOL) as i32;
                        c.tries = ((f64::from(nrt) + f64::from(nc1)) * f64::from(nct + others)
                            / f64::from(nn1 + others + nc1)
                            - f64::from(c.value)
                            + 1.) as i32;
                        c.row_left = nrt - c.value;
                    }
                    // The last cell takes the rest of the row.
                    let last = cells[j].row_left;
                    prefix[ncols] = prefix[j + 1] + fact[last];
                    cells[ncols - 1].value = last;

                    let v = base + prefix[ncols];
                    if rows_left == 2 {
                        // One row left: it is determined.
                        let left = |i: usize| fact[cols[i] - cells[i].value];
                        let v = (2..ncols).fold(v + (left(0) + left(1)), |v, i| v + left(i));
                        best = c_fmin(v, best);
                    } else if rows_left == 3 && ncols == 2 {
                        // Two rows of two columns left: the 2 x 2 mode.
                        let nn1 = total - rows[first] + 2;
                        let ic1 = cols[0] - cells[0].value;
                        let ic2 = cols[1] - cells[1].value;
                        let r = rows[first + 1];
                        let n11 = (r + 1).wrapping_mul(ic1 + 1) / nn1;
                        let n12 = r - n11;
                        let v = v + (fact[n11] + fact[n12] + fact[ic1 - n11] + fact[ic2 - n12]);
                        best = c_fmin(v, best);
                    } else {
                        // The remaining column totals form a node of the next row.
                        let scratch = &mut *key;
                        scratch.clear();
                        scratch.extend(
                            cols.iter()
                                .zip(cells.iter())
                                .map(|(&c, cell)| (c - cell.value).max(0)),
                        );
                        sort_ascending(scratch)?;
                        let key = scratch
                            .iter()
                            .try_fold(0, |key: Key, &x| {
                                key.checked_mul(radix)?.checked_add(Key::from(x))
                            })
                            .ok_or(Error::SearchKeyOverflow)?;
                        stack.push(key, v)?;
                    }
                }
                Step::PopNode => {
                    if let Some((mut key, value)) = stack.pop() {
                        base = value;
                        for c in cols[1..].iter_mut().rev() {
                            *c = (key % radix) as i32;
                            key /= radix;
                        }
                        cols[0] = key as i32;
                        set_cols_after(cells, cols, total);

                        // Closed form if possible. (The quick test compares
                        // with the original's `iro[nro]` rather than the last
                        // remaining row.)
                        let remaining = &rows[first..first + rows_left];
                        let quick = rows[rows_left - 1] <= rows[first] + ncols as i32;
                        match closed_form_with(remaining, cols, quick, base, fact, scratch) {
                            Some(v) => {
                                base = v;
                                best = c_fmin(v, best);
                            }
                            None => step = Step::NewNode,
                        }
                    } else if rows_left > 2 && stack.pushed > 0 {
                        // Row `first` done: continue with the next one.
                        stack.next_level();
                        total -= rows[first];
                        first += 1;
                        rows_left -= 1;
                    } else {
                        return Ok(-best);
                    }
                }
            }
        }
    }
}

/// `cells[j].cols_after = total - cols[0] - ... - cols[j]`.
fn set_cols_after(cells: &mut [Cell], cols: &[i32], total: i32) {
    let mut left = total;
    for (cell, &c) in cells.iter_mut().zip(cols) {
        left -= c;
        cell.cols_after = left;
    }
}

/// Tables with at most one row or column, or 2 x 2.
fn trivial(irow: &[i32], icol: &[i32], ntot: i32, fact: &LogFactorials) -> Option<f64> {
    let neg_sum = |xs: &[i32]| xs.iter().fold(0., |lp, &x| lp - fact[x]);
    match (irow.len(), icol.len()) {
        (0, _) | (_, 0) => Some(0.),
        (1, _) => Some(neg_sum(icol)),
        (_, 1) => Some(neg_sum(irow)),
        (2, 2) => {
            // The mode of the hypergeometric distribution.
            let n11 = (irow[0] + 1).wrapping_mul(icol[0] + 1) / (ntot + 2);
            let n12 = irow[0] - n11;
            Some(-(fact[n11] + fact[n12] + fact[icol[0] - n11] + fact[icol[1] - n12]))
        }
        _ => None,
    }
}

/// [`special_table`] for either orientation, each tried only if the
/// orientation's quick test passes (the row test is passed in, as the
/// caller's differs from the obvious one in one place).
fn closed_form_with(
    rows: &[i32],
    cols: &[i32],
    quick_rows: bool,
    start: f64,
    fact: &LogFactorials,
    scratch: &mut SpecialScratch,
) -> Option<f64> {
    let quick_cols = cols[cols.len() - 1] <= cols[0] + rows.len() as i32;
    quick_rows
        .then(|| special_table(rows, cols, start, fact, scratch))
        .flatten()
        .or_else(|| {
            quick_cols
                .then(|| special_table(cols, rows, start, fact, scratch))
                .flatten()
        })
}

fn closed_form(
    rows: &[i32],
    cols: &[i32],
    start: f64,
    fact: &LogFactorials,
    scratch: &mut SpecialScratch,
) -> Option<f64> {
    let quick_rows = rows[rows.len() - 1] <= rows[0] + cols.len() as i32;
    closed_form_with(rows, cols, quick_rows, start, fact, scratch)
}

/// Closed-form optimum for "special" tables (C: `f10act`), where each
/// column can be spread as evenly as possible over the rows: column `c`
/// gets `c / nrow` in every row plus one more in `c % nrow` of them. If the
/// row totals allow that, returns `start` plus the `sum ln(cell!)` of it.
fn special_table(
    rows: &[i32],
    cols: &[i32],
    start: f64,
    fact: &LogFactorials,
    scratch: &mut SpecialScratch,
) -> Option<f64> {
    let SpecialScratch { extra, split } = scratch;
    let nrow = rows.len() as i32;
    split.clear();
    split.extend(cols.iter().map(|&c| (c / nrow, c % nrow)));
    // extra[i]: number of columns giving an extra cell to more than i rows.
    extra.clear();
    extra.resize(rows.len() - 1, 0);
    let mut even = 0; // cells every row gets
    for &(q, m) in split.iter() {
        even += q;
        if m != 0 {
            extra[(m - 1) as usize] += 1;
        }
    }
    for i in (0..extra.len().saturating_sub(1)).rev() {
        extra[i] += extra[i + 1];
    }

    // Feasible iff, for every k, the k largest rows can absorb their share.
    let mut slack = 0;
    for (&e, &r) in extra.iter().zip(rows[1..].iter().rev()) {
        slack += even + e - r;
        if slack < 0 {
            return None;
        }
    }

    Some(split.iter().fold(start, |v, &(q, m)| {
        v + (f64::from(m) * fact[q + 1] + f64::from(nrow - m) * fact[q])
    }))
}

/// Reusable buffers for [`special_table`]: per-row extra-cell counts and
/// each column's (quotient, remainder) over the rows.
#[derive(Default)]
struct SpecialScratch {
    extra: Vec<i32>,
    split: Vec<(i32, i32)>,
}
