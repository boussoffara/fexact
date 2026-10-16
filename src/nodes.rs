//! Enumerating the daughters of a node (C: `f7xact` and the start of
//! `f2xact`'s outer loop).
//!
//! A node is a vector of row totals `rows` (sorted ascending) still to be
//! filled. Going to a daughter fills the next column, whose total is `n`:
//! `column[i] <= rows[i]` cells go to row `i`, and the daughter's row totals
//! are `rows - column`. [`ColumnFill`] enumerates every such `column`, in
//! the order of the original.

/// Cursor over the ways of filling a column.
pub(crate) struct ColumnFill {
    /// The row the next cell is taken from (C: `kd`).
    k: usize,
    /// The lowest row that can take another cell, once known (C: `ks`).
    ks: Option<usize>,
    /// Only rows `..len` are ever filled (C: `kmax`).
    len: usize,
}

impl ColumnFill {
    /// The first filling: greedily put the column total into the largest
    /// rows. Returns `None` if the rows cannot hold `n` cells.
    pub(crate) fn first(n: i32, rows: &[i32], column: &mut [i32]) -> Option<Self> {
        column.fill(0);
        let mut n = n;
        let mut k = rows.len();
        let mut len = rows.len();
        loop {
            k -= 1;
            column[k] = n.min(rows[k]);
            if column[len - 1] == 0 {
                len -= 1;
            }
            n -= column[k];
            if n == 0 || k == 0 {
                break;
            }
        }
        (n == 0).then_some(ColumnFill { k, ks: None, len })
    }

    /// Advances `column` to the next filling. Returns false when there are
    /// no more.
    pub(crate) fn advance(&mut self, rows: &[i32], column: &mut [i32]) -> bool {
        let ks = *self.ks.get_or_insert_with(|| {
            (0..)
                .find(|&i| column[i] != rows[i])
                .expect("a row with room")
        });
        // (Locals rather than fields: this runs once per daughter.)
        let mut k = self.k;

        if column[k] > 0 && k > ks {
            // Move one cell from row k to the highest lower row with room.
            column[k] -= 1;
            k -= 1;
            while rows[k] == 0 {
                k -= 1;
            }
            let mut m = k;
            while column[m] >= rows[m] {
                m -= 1;
            }
            column[m] += 1;
            if m == ks && column[m] == rows[m] {
                self.ks = Some(k);
            }
            self.k = k;
            return true;
        }

        loop {
            // Take one cell from the next row above k that has any...
            let Some(kk) = (k + 1..self.len).find(|&i| column[i] > 0) else {
                return false;
            };
            // ...and refill rows below kk greedily from the top with it plus
            // everything at or below k.
            let mut mm = 1 + column[..=k].iter().sum::<i32>();
            column[..=k].fill(0);
            k = kk;
            loop {
                k -= 1;
                column[k] = mm.min(rows[k]);
                mm -= column[k];
                if mm == 0 || k == 0 {
                    break;
                }
            }
            if mm > 0 {
                // The rows below kk cannot hold them: try the next row up.
                if kk + 1 == self.len {
                    return false;
                }
                k = kk;
                continue;
            }
            column[kk] -= 1;
            self.ks = Some((0..=k).find(|&i| column[i] < rows[i]).unwrap_or(k + 1));
            self.k = k;
            return true;
        }
    }
}
