# Changelog

## 0.2.0 (2026-10-16)

`fexact` is now written in Rust and built with maturin. Wherever the C
extension was correct, `fexact` returns the same p-values to the bit and the
same errors (checked on ~15k recorded cases), about twice as fast.

### New

- `fexact_many(arrays, workspace=200000, hybrid=None, threads=None)`: p-values
  of many tables, computed in parallel without the GIL. Raises
  `ValueError("table i: ...")` for the first table that fails.
- `fexact_simulate(array, samples=10000, seed=None, threads=None)`: Monte Carlo
  p-value, like R's `fisher.test(simulate.p.value = TRUE)`. Multithreaded,
  reproducible for a given seed whatever the number of threads, with no
  workspace limit.
- A Rust API (`fisher_exact`, `fisher_exact_many`, `fisher_simulated`, and
  `compat::fexact` with the C signature).

### Fixed

- Wrong p-values on wide tables: the C code's 32-bit path counts overflowed
  silently, e.g. `[[10, 6, 7] * 7, [5, 7, 9] * 7]` gave 0.781 instead of
  0.874 (and 0.329 instead of 0.929 with 30 columns). They are now exact
  (checked against Monte Carlo simulation and R).
- Tables with large margins that failed with "FEXACT error 501" (hash key
  overflow), such as `[[100000, 6, 7], [100000, 7, 9]]`, are now computed:
  node keys are 64-bit.
- Counts beyond 2^31 were silently wrapped to 32 bits (`2**32 + 5` was read
  as `5`). They are now read as 64-bit: `fexact` reports `Out of workspace.`
  for them, and `fexact_simulate` handles them.

### Changed

- Some tables that raised `ValueError` now return p-values (see Fixed).
- The error "bug negative key" is now "integer overflow in exact computation"
  (R's message), raised only if a search key would overflow 64 bits.
- Building from source needs a Rust toolchain; wheels are provided for Linux
  (x86_64, aarch64), macOS (x86_64, arm64) and Windows (x64), for Python 3.9
  and later.
- The license metadata is `GPL-3.0-only` (unchanged terms).

## 0.1.2

C extension (R's `fexact.c`).
