//! Nodes of the network and the past path lengths that reach them
//! (C: the `key`/`ipoin`/`stp`/`ifrq` tables, `f5xact` "PUT" and `f6xact`
//! "GET").
//!
//! Two stages are live at a time: the one whose nodes are being expanded
//! (read) and the next one, collecting their daughters (write). Each node is
//! a slot of an open-addressing hash table keyed by its encoded row totals.
//! The distinct past path lengths reaching a node (merged when within `tol`)
//! are entries of a per-stage pool, linked into a list for iteration and into
//! a binary search tree for lookup. Insertion order and tree shape follow the
//! C code exactly, since both affect the result.

use crate::error::Error;
use crate::{slot_of, Key};

/// Marks an unused node slot.
const EMPTY: Key = -1;
/// Marks a missing entry link.
const NIL: u32 = u32::MAX;

fn link(id: u32) -> Option<usize> {
    (id != NIL).then_some(id as usize)
}

/// Past path lengths, stored column-wise; grows as paths are added.
#[derive(Default)]
struct Pool {
    lengths: Vec<f64>,
    /// Number of paths (C: `int`, which overflows on wide tables; f64 is
    /// exact up to 2^53).
    freqs: Vec<f64>,
    /// Next entry of the same node.
    next: Vec<u32>,
    /// Search-tree children.
    left: Vec<u32>,
    right: Vec<u32>,
}

impl Pool {
    fn push(&mut self, length: f64, freq: f64, next: u32) -> u32 {
        let id = self.lengths.len() as u32;
        self.lengths.push(length);
        self.freqs.push(freq);
        self.next.push(next);
        self.left.push(NIL);
        self.right.push(NIL);
        id
    }

    fn clear(&mut self) {
        self.lengths.clear();
        self.freqs.clear();
        self.next.clear();
        self.left.clear();
        self.right.clear();
    }
}

/// One stage of the network: its nodes and their past paths.
pub(crate) struct Stage {
    /// Encoded row totals of the node in each slot, or [`EMPTY`].
    keys: Vec<Key>,
    /// First (root) entry of each slot's past paths.
    heads: Vec<u32>,
    pool: Pool,
    /// Maximum number of entries (C: `ldstp`).
    capacity: usize,
    /// Slot of the node targeted by the last [`put`](Self::put) (C's
    /// `static itp` in `f5xact`).
    last_slot: usize,
}

impl Stage {
    fn new(slots: usize, capacity: usize) -> Self {
        Stage {
            keys: vec![EMPTY; slots],
            heads: vec![NIL; slots],
            pool: Pool::default(),
            capacity,
            last_slot: 0,
        }
    }

    /// The (length, frequency) pairs of the past paths of the node in `slot`.
    pub(crate) fn paths(&self, slot: usize) -> impl Iterator<Item = (f64, f64)> + '_ {
        let pool = &self.pool;
        std::iter::successors(link(self.heads[slot]), |&i| link(pool.next[i]))
            .map(|i| (pool.lengths[i], pool.freqs[i]))
    }

    fn add_entry(&mut self, length: f64, freq: f64, next: u32) -> Result<u32, Error> {
        if self.pool.lengths.len() == self.capacity {
            return Err(Error::PastPathTableFull);
        }
        Ok(self.pool.push(length, freq, next))
    }

    /// Records `freq` paths of length `pastp` reaching node `key`. With
    /// `lookup == false` the node is the one of the previous call (saves a
    /// hash lookup when putting several paths of one node).
    pub(crate) fn put(
        &mut self,
        key: Key,
        pastp: f64,
        freq: f64,
        lookup: bool,
        tol: f64,
    ) -> Result<(), Error> {
        if lookup {
            // Linear probing from key % slots, wrapping around.
            let n = self.keys.len();
            let start = slot_of(key, n);
            let slot = (start..n)
                .chain(0..start)
                .find(|&s| self.keys[s] == key || self.keys[s] == EMPTY)
                .ok_or(Error::NodeTableFull)?;
            self.last_slot = slot;

            if self.keys[slot] == EMPTY {
                // New node: its first entry is the root of its tree.
                self.keys[slot] = key;
                self.heads[slot] = self.add_entry(pastp, freq, NIL)?;
                return Ok(());
            }
        }

        // Look for a length within `tol` of pastp in the node's tree.
        let (lo, hi) = (pastp - tol, pastp + tol);
        let pool = &mut self.pool;
        let mut i = self.heads[self.last_slot] as usize;
        let went_left = loop {
            let (child, left) = if pool.lengths[i] < lo {
                (pool.left[i], true)
            } else if pool.lengths[i] > hi {
                (pool.right[i], false)
            } else {
                pool.freqs[i] += freq;
                return Ok(());
            };
            if child == NIL {
                break left;
            }
            i = child as usize;
        };

        // Not found: add an entry as a child of the last node visited (`i`),
        // and after it in the node's list.
        let id = self.add_entry(pastp, freq, self.pool.next[i])?;
        self.pool.next[i] = id;
        if went_left {
            self.pool.left[i] = id;
        } else {
            self.pool.right[i] = id;
        }
        Ok(())
    }
}

pub(crate) struct PastPaths {
    stages: [Stage; 2],
    /// Index into `stages` of the read stage.
    read: usize,
    /// Next read-stage slot [`pop`](Self::pop) examines (C: `last`).
    cursor: usize,
}

impl PastPaths {
    /// Tables with `slots` node slots and room for `capacity` past paths per
    /// stage, holding just the root node (one path of length 0) in slot 0 of
    /// the read stage. The root is not in the key table, so it is never
    /// popped: the caller starts with it.
    pub(crate) fn new(slots: usize, capacity: usize) -> Self {
        let mut root = Stage::new(slots, capacity);
        root.heads[0] = root.pool.push(0., 1., NIL);
        PastPaths {
            stages: [root, Stage::new(slots, capacity)],
            read: 0,
            cursor: 0,
        }
    }

    /// The read stage and the write stage.
    pub(crate) fn split(&mut self) -> (&Stage, &mut Stage) {
        let [a, b] = &mut self.stages;
        if self.read == 0 {
            (a, b)
        } else {
            (b, a)
        }
    }

    /// Removes the next node of the read stage, decoding its row totals into
    /// `rows` (key = sum of `rows[i] * radix[i]`), and returns its slot.
    /// Returns `None` once the stage is exhausted.
    pub(crate) fn pop(&mut self, rows: &mut [i32], radix: &[Key]) -> Option<usize> {
        let st = &mut self.stages[self.read];
        while self.cursor < st.keys.len() {
            let slot = self.cursor;
            self.cursor += 1;
            let mut key = std::mem::replace(&mut st.keys[slot], EMPTY);
            if key == EMPTY {
                continue;
            }
            // Each digit is at most a row total, so fits an i32.
            for (r, &m) in rows.iter_mut().zip(radix).skip(1).rev() {
                *r = (key / m) as i32;
                key %= m;
            }
            rows[0] = key as i32;
            return Some(slot);
        }
        self.cursor = 0;
        None
    }

    /// Makes the write stage the read stage, with an empty write stage.
    /// (Its slots were all emptied by `pop`.)
    pub(crate) fn next_stage(&mut self) {
        self.stages[self.read].pool.clear();
        self.read = 1 - self.read;
    }
}
