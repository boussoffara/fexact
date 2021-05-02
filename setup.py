#!/usr/bin/env python

from distutils.core import setup, Extension
import numpy
setup(name             = "fexact",
      version          = "1.0",
      description      = "Example code for blog post.",
      author           = "J. David Lee",
      author_email     = "contact@crumpington.com",
      maintainer       = "contact@crumpington.com",
      url              = "https://www.crumpington.com",
      ext_modules      = [
          Extension(
              'fexact', ['fisher_acm/fexact_wrapper.c','fisher_acm/fexact.c','fisher_acm/helpers.c'],include_dirs=[numpy.get_include()],
              extra_compile_args=[]),

      ],

)
