# How to contribute

## Setup

`fexact` is a Rust crate (`src/`) with Python bindings built by
[maturin](https://www.maturin.rs/). You need a Rust toolchain
(<https://rustup.rs>) and Python >= 3.9.

```bash
python -m venv .venv && source .venv/bin/activate
pip install maturin pytest numpy
maturin develop --release    # builds and installs the `fexact` module
```

## Tests

```bash
cargo test     # Rust unit tests + tests/golden.rs (~15k cases recorded from the C code)
pytest         # Python tests
```

`tests/test_parity.py` compares against the original C extension, kept in
`legacy/`. It is skipped unless that is built:

```bash
pip install setuptools
python legacy/setup.py build_ext --inplace
```

The exact test must stay bit-identical to the C code wherever the C code is
correct: `tests/golden.rs` checks this. If you change behaviour on purpose,
explain why in the pull request and in `CHANGELOG.md`.

## Codestyle

```bash
cargo fmt
cargo clippy --features python --all-targets -- -D warnings
```

CI runs both, along with the tests on Linux, macOS and Windows.

## Before submitting

1. Add tests for your changes.
1. Run `cargo fmt`, clippy, `cargo test` and `pytest`.
1. Update `README.md` and `CHANGELOG.md` if you changed something user-visible.

## Other help

You can contribute by spreading a word about this library.
It would also be a huge contribution to write
a short article on how you are using this project.
You can also share your best practices with us.
