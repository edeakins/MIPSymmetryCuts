from setuptools import setup
from Cython.Build import cythonize
from distutils.extension import Extension
extensions = [
    Extension("inner_build_lp_loop", ["inner_build_lp_loop.pyx"]
              , language="c++"
              )
]

setup (
    ext_modules = cythonize (extensions)
)