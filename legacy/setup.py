#!/usr/bin/env python
"""Builds the original C extension, kept as the reference for the Rust port.

    python legacy/setup.py build_ext --inplace

puts `fexact.*.so` next to this file, where `tests/c_reference.py` finds it
(`tests/test_parity.py`, `tests/gen_golden.py`). It is not installed: the
`fexact` package is built from the Rust crate with maturin.
"""
import os
import sys

import numpy
from setuptools import Extension, setup

here = os.path.dirname(os.path.abspath(__file__))
os.chdir(here)

# fexact_wrapper.c passes an incompatible pointer type, an error since GCC 14.
extra = [] if sys.platform == "win32" else ["-Wno-error=incompatible-pointer-types"]

setup(
    name="fexact-c-reference",
    ext_modules=[
        Extension(
            "fexact",
            ["fexact_wrapper.c", "fexact.c", "helpers.c"],
            include_dirs=[numpy.get_include(), "."],
            extra_compile_args=extra,
        ),
    ],
)
