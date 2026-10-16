//! The network algorithm (C: `f2xact`).
//!
//! Tables with the observed margins are built column by column. A node at
//! stage `k` is the (sorted) vector of row totals left for the last `k`
//! columns; an arc fills one column, and the arc lengths along a path sum to
//! the log probability of the table it builds (up to constants). For every
//! node the algorithm keeps the distinct lengths of the paths reaching it
//! ("past paths", see [`PastPaths`]) with their multiplicities.
//!
//! Expanding a node, each past path is classified using bounds on the
//! remaining path length ([`LongestPath`], [`ShortestPath`]):
//! - every completion is at most as probable as the observed table: all of
//!   them count towards the p-value at once;
//! - no completion is: the path is dropped;
//! - otherwise it is pushed to the daughter node (or, in hybrid mode, its
//!   completions are approximated by a chi-squared tail when the daughter's
//!   expected counts satisfy the Cochran condition).
//!
//! Path lengths are kept with the C code's exact sequence of floating-point
//! operations: the goal is bit-identical p-values.

use crate::error::Error;
use crate::longest_path::LongestPath;
use crate::math::{c_fmax, gamma_cdf, log_multinomial, LogFactorials};
use crate::nodes::ColumnFill;
use crate::past_paths::{PastPaths, Stage};
use crate::shortest_path::ShortestPath;
use crate::sort::sort_ascending;
use crate::workspace::Layout;
use crate::{slot_of, Hybrid, Key, Output, Table, TOL};

/// `ln(2 pi)`, as written in the C code.
#[allow(clippy::excessive_precision)]
const LN_2PI: f64 = 1.83787706640934548356065947281;

pub(crate) fn network(
    table: &Table,
    hybrid: Option<&Hybrid>,
    layout: &Layout,
) -> Result<Output, Error> {
    if table.ncols() <= 1 {
        return Err(Error::TooFewColumns);
    }
    let (mut rows, cols, transposed) = margins(table)?;
    let radix = key_radix(&rows)?;
    let ntot: i32 = rows.iter().sum();
    let fact = LogFactorials::new(ntot);
    let obs = observed_length(table, rows.len(), &cols, transposed, &fact);
    let dro = log_multinomial(ntot, &rows, &fact);

    let nco = cols.len();
    let mut search = Search {
        obs,
        dro,
        hybrid,
        bounds: BoundsCache::new(2 * layout.ldkey),
        longest: LongestPath::new(nco, layout.ldst),
        shortest: ShortestPath::new(table.nrows() + table.ncols() + 1, nco),
        pending: None,
        pre: 0.,
        daughter: vec![0; rows.len()],
        fact,
        radix,
    };
    let mut paths = PastPaths::new(layout.ldkey, layout.ldstp);

    // Stage k fills column nco - k; the last column is implied.
    for k in (2..=nco).rev() {
        let cols = &cols[nco - k..];
        // The first stage has just the root, in slot 0.
        let mut node = if k == nco {
            Some(0)
        } else {
            paths.pop(&mut rows, &search.radix)
        };
        while let Some(slot) = node {
            let (read, write) = paths.split();
            search.expand(&rows, cols, read, slot, write)?;
            node = paths.pop(&mut rows, &search.radix);
        }
        paths.next_stage();
        search.bounds.clear();
    }

    match search.pending {
        Some(e) => Err(e),
        None => Ok(Output {
            p_value: search.pre,
            table_probability: (obs - dro).exp(),
        }),
    }
}

/// Sorted row and column totals, transposed if needed so that there are no
/// more rows than columns (the flag says whether they were).
fn margins(table: &Table) -> Result<(Vec<i32>, Vec<i32>, bool), Error> {
    let (nrow, ncol) = (table.nrows(), table.ncols());
    let mut rows: Vec<i32> = (0..nrow)
        .map(|i| (0..ncol).map(|j| table.cell(i, j)).sum())
        .collect();
    let mut cols: Vec<i32> = (0..ncol)
        .map(|j| (0..nrow).map(|i| table.cell(i, j)).sum())
        .collect();
    sort_ascending(&mut rows)?;
    sort_ascending(&mut cols)?;
    let transposed = nrow > ncol;
    if transposed {
        std::mem::swap(&mut rows, &mut cols);
    }
    Ok((rows, cols, transposed))
}

/// Length of the observed path (C: `obs`), offset by TOL so that ties with
/// the observed table count as "at least as extreme". The terms are added in
/// the C code's order: per column of the table as stored, paired with the
/// sorted column totals.
fn observed_length(
    table: &Table,
    nro: usize,
    cols: &[i32],
    transposed: bool,
    fact: &LogFactorials,
) -> f64 {
    let cell = |i, j| {
        if transposed {
            table.cell(j, i)
        } else {
            table.cell(i, j)
        }
    };
    cols.iter().enumerate().fold(TOL, |obs, (j, &c)| {
        let dd = (0..nro).fold(0., |dd, i| dd + fact[cell(i, j)]);
        obs + (fact[c] - dd)
    })
}

/// Multipliers encoding a row-total vector `r` (each `r[i] <= rows[i]`)
/// as the integer `sum r[i] * radix[i]`.
fn key_radix(rows: &[i32]) -> Result<Vec<Key>, Error> {
    let mut radix = Vec::with_capacity(rows.len());
    let mut m: Key = 1;
    for &r in rows {
        if Key::from(r) + 1 > Key::MAX / m {
            return Err(Error::KeyOverflow);
        }
        radix.push(m);
        m *= Key::from(r) + 1;
    }
    Ok(radix)
}

fn encode(rows: &[i32], radix: &[Key]) -> Key {
    rows.iter()
        .zip(radix)
        .map(|(&r, &m)| Key::from(r) * m)
        .sum()
}

/// State shared by all node expansions.
struct Search<'a> {
    fact: LogFactorials,
    radix: Vec<Key>,
    /// Observed path length, and the log multinomial of the row totals.
    obs: f64,
    dro: f64,
    hybrid: Option<&'a Hybrid>,
    bounds: BoundsCache,
    longest: LongestPath,
    shortest: ShortestPath,
    /// C records a full bounds cache but carries on until the next past
    /// path is stored, whose own error (if any) takes precedence.
    pending: Option<Error>,
    /// The p-value accumulated so far (C: `pre`).
    pre: f64,
    /// Scratch: the daughter's row totals.
    daughter: Vec<i32>,
}

/// How the past paths of a node are classified on the way to a daughter.
struct Arc {
    /// The daughter's node key.
    key: Key,
    /// Arc length (C: `ddf`), and that plus the daughter's multinomial term
    /// (C: `drn`).
    length: f64,
    drn: f64,
    /// Past paths up to `all_extreme` long have only extreme completions;
    /// those from `none_extreme` on have none (C: `obs3`, `obs2`).
    all_extreme: f64,
    none_extreme: f64,
    /// Chi-squared offset and degrees of freedom, in hybrid mode.
    chisq: Option<(f64, f64)>,
}

impl Search<'_> {
    /// Expands the node in `slot` of the `read` stage, whose row totals are
    /// `rows`, by the first of `cols`; daughters go to the `write` stage.
    fn expand(
        &mut self,
        rows: &[i32],
        cols: &[i32],
        read: &Stage,
        slot: usize,
        write: &mut Stage,
    ) -> Result<(), Error> {
        let (&n, rest) = cols.split_first().expect("at least two columns");
        let n_rest: i32 = rest.iter().sum();
        let mut column = vec![0; rows.len()];
        let Some(mut fill) = ColumnFill::first(n, rows, &mut column) else {
            return Ok(());
        };
        loop {
            let arc = self.arc(rows, &column, n, rest, n_rest)?;
            self.follow(&arc, read.paths(slot), write)?;
            if !fill.advance(rows, &mut column) {
                return Ok(());
            }
        }
    }

    /// The arc to the daughter obtained by filling the column with `column`.
    fn arc(
        &mut self,
        rows: &[i32],
        column: &[i32],
        n: i32,
        rest: &[i32],
        n_rest: i32,
    ) -> Result<Arc, Error> {
        let fact = &self.fact;
        for ((d, &r), &c) in self.daughter.iter_mut().zip(rows).zip(column) {
            *d = r - c;
        }
        // With more than one column left the daughter is a node: sort it and
        // skip its zero rows. (Otherwise the remaining column is determined;
        // the summation order of `drn` then follows the unsorted totals.)
        let nonzero = if rest.len() > 1 {
            sort_ascending(&mut self.daughter)?;
            self.daughter
                .iter()
                .position(|&x| x != 0)
                .unwrap_or(self.daughter.len())
        } else {
            0
        };
        let sub_rows = &self.daughter[nonzero..];
        let length = log_multinomial(n, column, fact);
        let drn = log_multinomial(n_rest, sub_rows, fact) - self.dro + length;

        if rest.len() == 1 {
            // One completion: extreme or not.
            let t = self.obs - drn - self.dro;
            return Ok(Arc {
                key: 0,
                length,
                drn,
                all_extreme: t,
                none_extreme: t,
                chisq: None,
            });
        }

        let key = encode(&self.daughter, &self.radix);
        let rem = rest[2..].iter().fold(
            self.obs - fact[rest[0]] - fact[rest[1]] - length,
            |rem, &x| rem - fact[x],
        );
        let bounds = match self.bounds.find(key) {
            Lookup::Found(b) => b,
            Lookup::Full(b) => {
                self.pending = Some(Error::BoundsCacheFull);
                b
            }
            Lookup::Vacant(slot) => {
                let dspt = self.obs - rem - length;
                let longest = self.longest.compute(sub_rows, rest, n_rest, fact)?;
                let shortest = self.shortest.compute(sub_rows, rest, dspt, fact);
                let b = PathBounds {
                    longest: at_most_zero(longest),
                    shortest: at_most_zero(shortest),
                    chisq: self.hybrid.and_then(|h| {
                        chi_square_offset(sub_rows, rest, n_rest, self.obs - self.dro, fact, h)
                    }),
                };
                self.bounds.insert(slot, key, b);
                b
            }
        };
        let df = ((sub_rows.len() - 1) * (rest.len() - 1)) as f64;
        Ok(Arc {
            key,
            length,
            drn,
            all_extreme: rem - bounds.longest,
            none_extreme: rem - bounds.shortest,
            chisq: bounds.chisq.map(|tmp| (tmp, df)),
        })
    }

    /// Classifies past paths `(length, frequency)` along `arc`.
    fn follow(
        &mut self,
        arc: &Arc,
        paths: impl Iterator<Item = (f64, f64)>,
        write: &mut Stage,
    ) -> Result<(), Error> {
        let mut lookup = true;
        for (pastp, freq) in paths {
            if pastp <= arc.all_extreme {
                self.pre += freq * (pastp + arc.drn).exp();
            } else if pastp < arc.none_extreme {
                if let Some((tmp, df)) = arc.chisq {
                    let x = c_fmax(0., tmp + (pastp + arc.drn) * 2.) / 2.;
                    let pv = 1. - gamma_cdf(x, df / 2.);
                    self.pre += freq * (pastp + arc.drn).exp() * pv;
                } else {
                    // Undecided: carry the path to the daughter.
                    write.put(arc.key, pastp + arc.length, freq, lookup, TOL)?;
                    if let Some(e) = self.pending {
                        return Err(e);
                    }
                    lookup = false;
                }
            }
        }
        Ok(())
    }
}

/// C's `if (x > 0.) x = 0.;` (unlike `f64::min`, keeps NaN and -0.0).
fn at_most_zero(x: f64) -> f64 {
    if x > 0. {
        0.
    } else {
        x
    }
}

/// Bounds on the remaining path length from a node, as offsets (see
/// [`network`]), plus the chi-squared offset if the hybrid approximation
/// applies to it.
#[derive(Debug, Clone, Copy)]
struct PathBounds {
    longest: f64,
    shortest: f64,
    chisq: Option<f64>,
}

/// Marks an unused cache slot.
const NIL: u32 = u32::MAX;

/// Per-stage cache of [`PathBounds`] by node key, with open addressing
/// (C: `key2`, `LP`, `SP`, `tm`). Entries are stored densely so that
/// clearing costs only as much as the stage used.
struct BoundsCache {
    /// Index into `entries` for each slot, or [`NIL`].
    slots: Vec<u32>,
    /// (slot, key, bounds) in insertion order.
    entries: Vec<(usize, Key, PathBounds)>,
}

enum Lookup {
    Found(PathBounds),
    Vacant(usize),
    /// The table is full; C carries on with the entry where the search
    /// started.
    Full(PathBounds),
}

impl BoundsCache {
    fn new(slots: usize) -> Self {
        BoundsCache {
            slots: vec![NIL; slots],
            entries: Vec::new(),
        }
    }

    fn find(&self, key: Key) -> Lookup {
        let n = self.slots.len();
        let start = slot_of(key, n);
        for s in (start..n).chain(0..start) {
            match self.slots[s] {
                NIL => return Lookup::Vacant(s),
                e if self.entries[e as usize].1 == key => {
                    return Lookup::Found(self.entries[e as usize].2)
                }
                _ => {}
            }
        }
        Lookup::Full(self.entries[self.slots[start] as usize].2)
    }

    fn insert(&mut self, slot: usize, key: Key, bounds: PathBounds) {
        self.slots[slot] = self.entries.len() as u32;
        self.entries.push((slot, key, bounds));
    }

    fn clear(&mut self) {
        for (slot, ..) in self.entries.drain(..) {
            self.slots[slot] = NIL;
        }
    }
}

/// The hybrid approximation's offset for a subtable with row totals `rows`
/// and column totals `cols` (summing to `n`), if its expected counts meet
/// the [`Hybrid`] criteria. `obs_dro` is the observed path length minus the
/// log multinomial of the row totals.
#[allow(clippy::neg_cmp_op_on_partial_ord)] // C's NaN behaviour
fn chi_square_offset(
    rows: &[i32],
    cols: &[i32],
    n: i32,
    obs_dro: f64,
    fact: &LogFactorials,
    hybrid: &Hybrid,
) -> Option<f64> {
    // `rows[0] * cols[0]` is the smallest expected count (times n); the
    // products wrap on overflow as in C.
    let expected = |r: i32, c: i32| f64::from(r.wrapping_mul(c));
    if !(expected(rows[0], cols[0]) > f64::from(n) * hybrid.emin) {
        return None;
    }
    let large = rows
        .iter()
        .flat_map(|&r| cols.iter().map(move |&c| (r, c)))
        .filter(|&(r, c)| expected(r, c) >= f64::from(n) * hybrid.expect)
        .count() as i32;
    let (nr, nc) = (rows.len() as i32, cols.len() as i32);
    if !(f64::from(large * 100) >= f64::from(nc * nr) * hybrid.percnt) {
        return None;
    }

    let mut tmp = 0.;
    for &r in rows {
        tmp += fact[r] - fact[r - 1];
    }
    tmp *= f64::from(nc - 1);
    for &c in cols {
        tmp += f64::from(nr - 1) * (fact[c] - fact[c - 1]);
    }
    let df = f64::from((nr - 1) * (nc - 1));
    tmp += df * LN_2PI;
    tmp -= f64::from(nr * nc - 1) * (fact[n] - fact[n - 1]);
    // C stores "not applicable" as the value -9876.
    Some(obs_dro * -2. - tmp).filter(|&t| t != -9876.)
}
