//! Workspace budget.
//!
//! The C code carves every array out of one buffer of `workspace` ints. Its
//! budget decides whether a problem fits ("Out of workspace.") and how large
//! the hash tables are; both are reproduced here so results and errors match.
//! The arrays themselves are allocated at their natural sizes elsewhere.

use crate::error::Error;

/// Hash table sizes derived from the workspace budget.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Layout {
    /// Slots per stage in the node and path-bound tables (C: `ldkey`).
    pub ldkey: usize,
    /// Past-path-length entries per stage (C: `ldstp`, `mult * ldkey`).
    pub ldstp: usize,
    /// Half the size of the longest-path search stack (C: `n2_stack`).
    pub ldst: usize,
}

impl Layout {
    /// Replays the C code's `iwork()` allocations for a table with `nrow` x
    /// `ncol` cells summing to `ntot`.
    pub(crate) fn new(
        workspace: i64,
        mult: i64,
        nrow: usize,
        ncol: usize,
        ntot: i32,
    ) -> Result<Self, Error> {
        let iwkmax = 2 * (workspace / 2);
        let ldst = 200.max(iwkmax / 1000);
        let (nco, nro) = (nrow.max(ncol) as i64, nrow.min(ncol) as i64);
        let k = (nrow + ncol + 1) as i64;

        // Sizes in ints; a double counts as two. (C rounds a double's start
        // index up to even but does not add that padding to the running
        // total, hence the `- 1` slack in `ldkey` below.)
        let fixed = 2 * (i64::from(ntot) + 1) // log factorials
            + 3 * nco + 2 * nro // marginals, key radices, column and daughter rows
            + (5 * k + 2 * k * nco).max(7 * nco + 4 * ldst) // int scratch
            + 2 * (nco + 1 + 2 * ldst).max(k); // double scratch

        // The rest is split into hash tables. Per unit of `ldkey`: 2 key and
        // 2 node-pointer ints, 2*mult doubles + 6*mult ints of past paths,
        // and 2 keys + 3*2 doubles of path bounds = 18 + 10*mult ints.
        let ldkey = (iwkmax - fixed) / (18 + 10 * mult) - 1;
        if mult * ldkey > i64::from(i32::MAX) {
            return Err(Error::IntegerOverflow);
        }
        let ldstp = mult * ldkey;
        let tables =
            2 * ldkey + 2 * ldkey + 2 * (2 * ldstp) + 6 * ldstp + 3 * 2 * (2 * ldkey) + 2 * ldkey;

        // C fails if any allocation overran the buffer. The hash tables are
        // sized to fit, so this reduces to the fixed part and the total. With
        // `ldkey <= 0` the C code divides by zero or corrupts its tables.
        if fixed > iwkmax || fixed + tables > iwkmax || ldkey <= 0 {
            return Err(Error::OutOfWorkspace);
        }
        Ok(Layout {
            ldkey: ldkey as usize,
            ldstp: ldstp as usize,
            ldst: ldst as usize,
        })
    }
}
