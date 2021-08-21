from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "misc_cython: equations for thermophysical proptery estimation",
    ext_modules = cythonize('misc_cython.pyx'),
)
