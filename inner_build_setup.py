from setuptools import setup
from Cython.Build import cythonize

setup (
    ext_modules = cythonize ("inner_build_lp_loop.pyx")
)