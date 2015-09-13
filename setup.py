"""
Project website: https://bitbucket.org/yymao/abundancematching
Copyright (c) 2015 Yao-Yuan Mao (yymao)
"""

from setuptools import setup, find_packages
import os
from subprocess import check_call

def make():
    here = os.path.abspath(os.path.dirname(__file__))
    cwd = os.getcwd()
    os.chdir(os.path.join(here, 'fiducial_deconvolute'))
    check_call(['make'])
    os.chdir(cwd)

make()

setup(
    name='AbundanceMatching',
    version='0.1.0',
    description='A python module to create (interpolate and extrapolate) abundance functions and also provide fiducial deconvolution (with Peter Behroozi\'s implmentation) and abundance matching.',
    url='https://bitbucket.org/yymao/abundancematching',
    author='Yao-Yuan Mao, Peter Behroozi',
    author_email='Yao-Yuan Mao <yymao@alumni.stanford.edu>',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 2 :: Only',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    use_2to3=True,
    packages=find_packages(),
    package_data={
        'fiducial_deconvolute': ['fiducial_deconvolute.so'],
    },
    install_requires = ['numpy','scipy'],
)
