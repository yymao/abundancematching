"""
Project website: https://bitbucket.org/yymao/abundancematching
Copyright (c) 2015--2017 Yao-Yuan Mao (yymao)
"""

from setuptools import setup
import os
from subprocess import check_call


_make_pre = 'gcc -D_BSD_SOURCE -D_POSIX_SOURCE -D_POSIX_C_SOURCE=200809L -D_SVID_SOURCE -D_DARWIN_C_SOURCE -Wall -fno-math-errno -fPIC -shared'.split()
_make_post = '-lm -O3 -std=c99'.split()
here = os.path.abspath(os.path.dirname(__file__))
cpath = os.path.join(here, 'AbundanceMatching', 'fiducial_deconvolute')
check_call(_make_pre + [cpath+'.c', '-o', cpath+'.so'] + _make_post)

with open(os.path.join(here, 'AbundanceMatching', 'version.py')) as f:
    exec(f.read())


setup(
    name='AbundanceMatching',
    version=__version__,
    description='A python module to create (interpolate and extrapolate) abundance functions and also provide fiducial deconvolution (with Peter Behroozi\'s implmentation) and abundance matching.',
    url='https://bitbucket.org/yymao/abundancematching',
    download_url = 'https://bitbucket.org/yymao/abundancematching/get/v0.2.1.tar.gz',
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
    package_data={
        'AbundanceMatching': ['fiducial_deconvolute.c', 'fiducial_deconvolute.so'],
    },
    install_requires = ['numpy','scipy'],
)
