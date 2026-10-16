use std::fmt;

/// Why a p-value could not be computed.
///
/// [`Display`](fmt::Display) gives the legacy C extension's message verbatim
/// (R's for [`Error::SearchKeyOverflow`]), which is also what the Python
/// binding raises as `ValueError`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Error {
    /// The workspace size is negative (C: `malloc` failure).
    AllocationFailed,
    /// `nrow` exceeds the table's leading dimension.
    LeadingDimension,
    /// A table entry is negative.
    NegativeEntry,
    /// Every table entry is zero.
    AllZero,
    /// The table has fewer than two columns.
    TooFewColumns,
    /// The workspace is too small for the table.
    OutOfWorkspace,
    /// `mult * ldkey` overflows an `int`.
    IntegerOverflow,
    /// Node keys would overflow 64 bits (C errors 5 and 501, at 32 bits).
    KeyOverflow,
    /// The table caching path-length bounds of a stage is full (C error 6 in
    /// `f2xact`).
    BoundsCacheFull,
    /// The table of nodes of a stage is full (C error 6 in `f5xact`).
    NodeTableFull,
    /// The table of past path lengths is full (C error 7).
    PastPathTableFull,
    /// A longest-path search key would overflow 64 bits. (The C code's
    /// `int` keys wrap around instead, giving wrong p-values or "bug
    /// negative key".)
    SearchKeyOverflow,
    /// The longest-path search stack is full (C error 30).
    MarginStackFull,
    /// The sort routine's fixed-size stack overflowed (only for dimensions
    /// above ~1000).
    SortStackOverflow,
}

impl Error {
    /// The C extension's message for this error.
    pub fn message(self) -> &'static str {
        match self {
            Error::AllocationFailed => "Memory allocation failed",
            Error::LeadingDimension => "NROW must be less than or equal to LDTABL.",
            Error::NegativeEntry => "All elements of TABLE must be nonnegative.",
            Error::AllZero => {
                "All elements of TABLE are zero.\nPRT and PRE are set to missing values."
            }
            Error::TooFewColumns => "NCOL must be at least 2",
            Error::OutOfWorkspace => "Out of workspace.",
            Error::IntegerOverflow => "integer overflow",
            Error::KeyOverflow => {
                "The hash table key cannot be computed because the largest key\n\
                 is larger than the largest representable int.\n\
                 The algorithm cannot proceed.\n\
                 Reduce the workspace, consider using 'simulate.p.value=TRUE' or another algorithm."
            }
            Error::BoundsCacheFull => "LDKEY is too small.\n",
            Error::NodeTableFull => "LDKEY is too small for this problem.\n",
            Error::PastPathTableFull => "LDSTP is too small for this problem.\n",
            Error::SearchKeyOverflow => "integer overflow in exact computation",
            Error::MarginStackFull => "increase workspace",
            Error::SortStackOverflow => "This should never occur, please report!",
        }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.message())
    }
}

impl std::error::Error for Error {}
