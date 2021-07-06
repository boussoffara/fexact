#!/usr/bin/env python
from distutils.core import setup, Extension
import numpy

setup(
    name="fexact",
    version="0.0.8",
    description="Exact Fisher test for nxm contingency matrix based on fortran acm 643 and R",
    license="GPLv3",
    author="A. Boussoffara",
    author_email="boussoffara@me.com",
    maintainer="boussoffara@me.com",
    url="https://github.com/boussoffara/fexact",
    ext_modules=[
        Extension(
            "fexact",
            [
                "src/fexact_wrapper.c",
                "src/fexact.c",
                "src/helpers.c",
            ],
            include_dirs=[numpy.get_include(), "src/helpers.h", "src/fexact.h"],
            extra_compile_args=[],
        ),
    ],
)
