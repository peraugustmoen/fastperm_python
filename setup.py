from setuptools import setup, Extension
import numpy as np

def main():
    setup(
        license='MIT',
        description = 'The Fastperm package implements the algorithm proposed by Christensen (Inference for Bayesian Nonparametric Models with Binary Response Data via Permutation Counting, Bayesian Analysis, 1(1), 1-26,2023) for computing the permanent of a block rectangular matrix.',
        author = 'Dennis Christensen, Per August Jarval Moen',                  
        author_email = 'Dennis.Christensen@ffi.no, pamoen@math.uio.no',      
        url = 'https://github.com/peraugustmoen/fastperm',  
        download_url = 'https://github.com/peraugustmoen/fastperm/archive/refs/tags/v_01.tar.gz',
        name="fastperm",
        packages=["fastperm"],
        version="0.1",
        install_requires=["numpy"],
        include_dirs=[np.get_include()],
        ext_modules=[Extension("fastperm", ["get_log_permanent.c", "help_functions.c", "get_alphabetagamma.c",\
                "methods_sparse.c", "methods.c", "xxhash.c"])],
        headers = ["xxhash.h", "header.h"])



if __name__ == "__main__":
    main()
