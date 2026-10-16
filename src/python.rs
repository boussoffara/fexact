//! Python bindings, a drop-in replacement for the C extension in
//! `legacy/fexact_wrapper.c`, plus a batch function.

use std::num::NonZeroUsize;

use numpy::{PyArray2, PyArrayDescrMethods, PyArrayMethods, PyUntypedArray, PyUntypedArrayMethods};
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;

use crate::{
    fisher_exact, fisher_exact_many, fisher_simulated, Error, Hybrid, Options, Simulation, Table,
};

const ERR_ARRAY: &str =
    "array must be two-dimensional and cast of type np.int64, (np.array(x,dtype=np.int64)";

/// A table's counts in row-major order, with its shape.
struct Counts {
    cells: Vec<i64>,
    nrows: usize,
    ncols: usize,
}

impl Counts {
    fn extract(array: &Bound<'_, PyAny>) -> PyResult<Self> {
        let untyped = array.cast::<PyUntypedArray>().map_err(|_| {
            PyTypeError::new_err(format!(
                "fexact() argument 1 must be numpy.ndarray, not {}",
                array
                    .get_type()
                    .name()
                    .map_or_else(|_| "?".into(), |n| n.to_string())
            ))
        })?;
        if untyped.ndim() != 2
            || !untyped
                .dtype()
                .is_equiv_to(&numpy::dtype::<i64>(array.py()))
        {
            return Err(PyValueError::new_err(ERR_ARRAY));
        }
        let array = untyped.cast::<PyArray2<i64>>()?.readonly();
        let view = array.as_array();
        let (nrows, ncols) = view.dim();
        // (The C wrapper narrowed each value to a C int, silently wrapping
        // values beyond 2^31.)
        let cells = view.iter().copied().collect();
        Ok(Counts {
            cells,
            nrows,
            ncols,
        })
    }

    fn table(&self) -> Table<'_> {
        Table::from_row_major(&self.cells, self.nrows, self.ncols)
    }
}

/// Options with the C wrapper's semantics for `workspace` and `hybrid`.
fn options(workspace: i32, hybrid: Option<&Bound<'_, PyAny>>) -> PyResult<Options> {
    // C: malloc() of a negative size fails; -1 rounds down to 0.
    if 2 * (workspace / 2) < 0 {
        return Err(value_error(Error::AllocationFailed));
    }
    let hybrid = match hybrid {
        Some(h) => h.is_truthy()?,
        None => false,
    };
    Ok(Options {
        workspace: workspace.max(0) as u32,
        hybrid: hybrid.then_some(Hybrid::COCHRAN),
        ..Options::default()
    })
}

fn value_error(e: Error) -> PyErr {
    PyValueError::new_err(e.to_string())
}

/// Fisher's exact test p-value for a 2-D ``np.int64`` contingency table.
///
/// ``hybrid=True`` uses the Cochran-condition hybrid approximation
/// (expect=5, percnt=80, emin=1). Raises ``ValueError`` on invalid input or
/// when ``workspace`` is too small.
#[pyfunction]
#[pyo3(signature = (array, workspace = 200_000, hybrid = None))]
fn fexact(
    py: Python<'_>,
    array: &Bound<'_, PyAny>,
    workspace: i32,
    hybrid: Option<&Bound<'_, PyAny>>,
) -> PyResult<f64> {
    let options = options(workspace, hybrid)?;
    let counts = Counts::extract(array)?;
    py.detach(|| fisher_exact(&counts.table(), &options))
        .map(|r| r.p_value)
        .map_err(value_error)
}

/// ``fexact`` for each table of ``arrays``, computed in parallel on
/// ``threads`` threads (default: all cores). Returns the p-values in order;
/// raises ``ValueError`` naming the first table that fails.
#[pyfunction]
#[pyo3(signature = (arrays, workspace = 200_000, hybrid = None, threads = None))]
fn fexact_many(
    py: Python<'_>,
    arrays: &Bound<'_, PyAny>,
    workspace: i32,
    hybrid: Option<&Bound<'_, PyAny>>,
    threads: Option<NonZeroUsize>,
) -> PyResult<Vec<f64>> {
    let options = options(workspace, hybrid)?;
    let counts = arrays
        .try_iter()?
        .map(|a| Counts::extract(&a?))
        .collect::<PyResult<Vec<_>>>()?;
    let results = py.detach(|| {
        let tables: Vec<Table> = counts.iter().map(Counts::table).collect();
        fisher_exact_many(&tables, &options, threads)
    });
    results
        .into_iter()
        .enumerate()
        .map(|(i, r)| {
            r.map(|o| o.p_value)
                .map_err(|e| PyValueError::new_err(format!("table {i}: {e}")))
        })
        .collect()
}

/// Monte Carlo estimate of the ``fexact`` p-value from ``samples`` random
/// tables with the margins of ``array`` (as R's ``simulate.p.value``):
/// ``(1 + hits) / (samples + 1)``, where hits are tables at most as probable
/// as the observed one. Runs on ``threads`` threads (default: all cores); the
/// result depends on ``seed`` (default: random) but not on ``threads``.
#[pyfunction]
#[pyo3(signature = (array, samples = 10_000, seed = None, threads = None))]
fn fexact_simulate(
    py: Python<'_>,
    array: &Bound<'_, PyAny>,
    samples: u64,
    seed: Option<u64>,
    threads: Option<NonZeroUsize>,
) -> PyResult<f64> {
    let counts = Counts::extract(array)?;
    let simulation = Simulation {
        samples,
        seed,
        threads,
    };
    py.detach(|| fisher_simulated(&counts.table(), &simulation))
        .map(|r| r.p_value)
        .map_err(value_error)
}

#[pymodule]
#[pyo3(name = "fexact")]
fn fexact_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(fexact, m)?)?;
    m.add_function(wrap_pyfunction!(fexact_many, m)?)?;
    m.add_function(wrap_pyfunction!(fexact_simulate, m)?)?;
    Ok(())
}
