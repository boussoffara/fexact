#!/usr/bin/env python

import os
from distutils.core import setup, Extension
import numpy

setup(name             = "fexact",
      version          = "0.0.3",
      description      = "Exact Fisher test for nxm contingency matrix based on fortran acm 643 and R",
      license='GPLv3',
      author           = "A. Boussoffara",
      author_email     = "boussoffara@me.com",
      maintainer       = "boussoffara@me.com",
      url              = "https://www.github.com/boussoffara/fisher-acm",
      ext_modules      = [
          Extension(
              'fexact', ['fisher_acm/fexact_wrapper.c','fisher_acm/fexact.c','fisher_acm/helpers.c',
                         ],include_dirs=[numpy.get_include(),'fisher_acm/helpers.h','fisher_acm/fexact.h'],
              extra_compile_args=[]),

      ],

)
