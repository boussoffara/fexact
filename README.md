# fexact

<div align="center">

[![Build status](https://github.com/boussoffara/fexact/actions/workflows/test.yml/badge.svg?branch=master&event=push)](https://github.com/boussoffara/fexact/actions/workflows/test.yml)
[![Python Version](https://img.shields.io/pypi/pyversions/fexact.svg)](https://pypi.org/project/fexact/)
[![PyPI](https://img.shields.io/pypi/v/fexact.svg)](https://pypi.org/project/fexact/)
[![Semantic Versions](https://img.shields.io/badge/%F0%9F%9A%80-semantic%20versions-informational.svg)](https://github.com/boussoffara/fexact/releases)
[![License](https://img.shields.io/github/license/boussoffara/fexact)](https://github.com/boussoffara/fexact/blob/master/LICENSE.txt)

`fexact` is Fisher's exact test for r x c contingency tables: a port of R's network algorithm (ACM 643), with a Monte Carlo variant

</div>

## Install

```bash
pip install fexact
```

Wheels are built for Linux, macOS and Windows (Python >= 3.9). Building from
source needs a Rust toolchain.

## Usage

```python
import numpy as np
from fexact import fexact, fexact_many, fexact_simulate

contingency = np.array([[5, 1, 5], [1, 1, 5]], dtype=np.int64)
fexact(contingency)                                  # 0.40045...
fexact(contingency, workspace=10**7, hybrid=True)    # Cochran-condition hybrid approximation

# Many tables in parallel (the GIL is released; threads=None uses all cores).
fexact_many([t1, t2, t3], workspace=10**6, threads=None)  # -> [p1, p2, p3]

# Monte Carlo p-value from random tables with the same margins, in parallel.
# Reproducible for a given seed, whatever the number of threads.
fexact_simulate(contingency, samples=100_000, seed=1, threads=None)
```

`fexact_simulate` works like R's `fisher.test(simulate.p.value = TRUE)`: it
returns `(1 + hits) / (samples + 1)`, where hits are the random tables at most
as probable as the observed one. It has no workspace limit, so it also handles
tables too large for the exact test, including counts beyond 2^31 (about 0.2 s
per million samples for a 10x10 table with n = 5000 on 32 threads).

Tables must be two-dimensional `np.int64` arrays. `fexact` raises `ValueError` on
invalid input or when `workspace` is too small; `fexact_many` raises
`ValueError("table i: ...")` for the first table that fails.

## Implementation

`fexact` is written in Rust (`src/`): a port of R's `fexact.c` (ACM TOMS 643,
Mehta & Patel's network algorithm) that returns bit-identical p-values and the
same errors as the original C extension, about twice as fast, except where the
C code's 32-bit integers overflow. There the C extension either fails (error
501) or, on wide tables such as 2 x 21 and above, silently returns wrong
p-values; this version uses 64-bit node keys and path counts and gets them
right (checked against R's and its own Monte Carlo simulation). `src/lib.rs` maps
each module to the C routine it ports. The C code is kept in `legacy/` as the
reference.

It can also be used as a Rust library:

```rust
use fexact::{fisher_exact, fisher_exact_many, Hybrid, Options, Table};

let counts = [5, 1, 5,
              1, 1, 5];
let table = Table::from_row_major(&counts, 2, 3);
let exact = fisher_exact(&table, &Options::default())?.p_value;
let hybrid = fisher_exact(&table, &Options { hybrid: Some(Hybrid::COCHRAN), ..Options::default() })?;
```

`fexact::fisher_simulated(&table, &Simulation { samples, seed, threads })`
is the Monte Carlo test. `fexact::compat::fexact` keeps the C signature (column-major table with leading
dimension, `expect`/`percnt`/`emin`, `workspace`, `mult`).

## Development

```bash
cargo test                       # unit tests + tests/golden.rs
pip install maturin pytest
maturin develop --release        # or: pip install .
pytest
```

- `tests/golden.rs` replays ~15k cases recorded from the C extension
  (`tests/data/golden.txt`) and requires exact agreement.
- `tests/test_parity.py` compares random tables against the C extension. It is
  skipped unless the C reference has been built with
  `python legacy/setup.py build_ext --inplace`, which is also what
  `tests/gen_golden.py` needs to regenerate the corpus.

## Contributing

See [CONTRIBUTING.md](https://github.com/boussoffara/fexact/blob/master/CONTRIBUTING.md).
Changes are listed in [CHANGELOG.md](https://github.com/boussoffara/fexact/blob/master/CHANGELOG.md).

## 📈 Releases

You can see the list of available releases on the [GitHub Releases](https://github.com/boussoffara/fexact/releases) page.

We follow [Semantic Versions](https://semver.org/) specification.

We use [`Release Drafter`](https://github.com/marketplace/actions/release-drafter). As pull requests are merged, a draft release is kept up-to-date listing the changes, ready to publish when you’re ready. With the categories option, you can categorize pull requests in release notes using labels.

For Pull Request this labels are configured, by default:

|               **Label**               |  **Title in Releases**  |
| :-----------------------------------: | :---------------------: |
|       `enhancement`, `feature`        |       🚀 Features       |
| `bug`, `refactoring`, `bugfix`, `fix` | 🔧 Fixes & Refactoring  |
|       `build`, `ci`, `testing`        | 📦 Build System & CI/CD |
|              `breaking`               |   💥 Breaking Changes   |
|            `documentation`            |    📝 Documentation     |
|            `dependencies`             | ⬆️ Dependencies updates |

You can update it in [`release-drafter.yml`](https://github.com/boussoffara/fexact/blob/master/.github/release-drafter.yml).

GitHub creates the `bug`, `enhancement`, and `documentation` labels for you. Dependabot creates the `dependencies` label. Create the remaining labels on the Issues tab of your GitHub repository, when you need them.

To release:

1. Bump the version with `bump2version patch|minor|major` (updates `Cargo.toml` and `Cargo.lock`, commits and tags), and push with `--follow-tags`.
2. Publish a pre-release on GitHub to upload to Test PyPI, then a release to upload to PyPI. Both build wheels for Linux, macOS and Windows plus an sdist.

## 🛡 License

[![License](https://img.shields.io/github/license/boussoffara/fexact)](https://github.com/boussoffara/fexact/blob/master/LICENSE.txt)

This project is licensed under the terms of the `GNU GPL v3.0` license (`GPL-3.0-only`). See [LICENSE.txt](https://github.com/boussoffara/fexact/blob/master/LICENSE.txt) for more details.

## 📃 Citation

```
@misc{fexact,
  author = {boussoffara},
  title = {`fexact`  is a port of R's Fisher test for nxm contingency tables (acm 643 algorithm) },
  year = {2021},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/boussoffara/fexact}}
}
```

## Credits

The algorithm is Mehta & Patel's network algorithm (ACM TOMS 643) as found in R's `fexact.c`. The project was started from [`python-package-template`](https://github.com/TezRomacH/python-package-template).
