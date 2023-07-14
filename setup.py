from setuptools import setup, Extension
import numpy
setup(
    name="fastperm",
    version="0.1",
    include_dirs=[numpy.get_include()],
    ext_modules=[Extension("fastperm", ["get_log_permanent.c", "help_functions.c", "get_alphabetagamma.c",\
                "methods_sparse.c", "methods.c", "xxhash.c"])]
)
