from setuptools import setup, Extension
import numpy as np

def main():
    setup(
        name="fastperm",
        version="0.1",
        include_dirs=[np.get_include()],
        ext_modules=[Extension("fastperm", ["log_sum_exp.c","get_log_permanent.c", "help_functions.c", "get_alphabetagamma.c",\
                "methods_sparse.c", "methods.c", "xxhash.c"])],
        headers = ["xxhash.h", "header.h"])

if __name__ == "__main__":
    main()
