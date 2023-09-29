from setuptools import setup, Extension
import numpy as np

def main():
    setup(
        #license='BSD-2-Clause',
        description = 'The perms package implements the algorithm proposed by Christensen (Inference for Bayesian Nonparametric Models with Binary Response Data via Permutation Counting, Bayesian Analysis, 1(1), 1-26,2023) for computing the permanent of a block rectangular matrix.',
        author = 'Dennis Christensen, Per August Jarval Moen',                  
        author_email = 'Dennis.Christensen@ffi.no, pamoen@math.uio.no',      
        url = 'https://github.com/peraugustmoen/perms',  
        download_url = 'https://github.com/peraugustmoen/fastperm/archive/refs/tags/v_1_12.tar.gz',
        name="fastperm",
        version="1.13",
        #install_requires=["numpy"],
        include_dirs=[np.get_include()],
        ext_modules=[Extension("perms", ["get_log_permanent.c", "help_functions.c", "get_alphabetagamma.c",\
                "methods_sparse.c", "methods.c", "xxhash.c"])],
        headers = ["xxhash.h", "header.h"])



if __name__ == "__main__":
    main()
