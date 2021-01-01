"""
Project website: https://github.com/yymao/abundancematching
Copyright (c) 2015-2020 Yao-Yuan Mao (yymao)
"""

import os
from setuptools import setup, Extension

module_fd = Extension(
    'AbundanceMatching._fiducial_deconvolute',
    libraries=['m'],
    extra_compile_args=["-std=c99", "-fno-math-errno"],
    sources=['AbundanceMatching/fiducial_deconvolute.c'],
)

with open(os.path.join(os.path.dirname(__file__), 'AbundanceMatching', 'version.py')) as f:
    exec(f.read())

setup(
    name='AbundanceMatching',
    version=__version__,
    description='A Python package for subhalo abundance matching (SHAM) with scatter. It creates abundance functions and runs Peter Behroozi\'s fiducial deconvolution code.',
    url='https://github.com/yymao/abundancematching',
    download_url='https://github.com/yymao/abundancematching/archive/v{}.tar.gz'.format(__version__),
    author='Yao-Yuan Mao',
    author_email='yymao.astro@gmail.com',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    use_2to3=False,
    packages=['AbundanceMatching'],
    install_requires=['numpy >=1.9.0', 'scipy >=0.15.0'],
    ext_modules=[module_fd],
)
