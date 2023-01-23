from setuptools import setup
from Cython.Build import cythonize
from distutils.extension import Extension
extensions = [
    Extension("inner_build_lp_loop", 
              sources=["inner_build_lp_loop.pyx"],
              extra_compile_args=['-O3'],
              language="c++"
              )
]

setup (
    ext_modules = cythonize (extensions)
)