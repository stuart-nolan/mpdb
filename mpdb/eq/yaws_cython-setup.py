from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "yaws_cython: yaws equations thermophysical proptery estimation",
    ext_modules = cythonize('yaws_cython.pyx'),
)
