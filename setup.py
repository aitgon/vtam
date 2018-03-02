import os
from setuptools import setup
from setuptools import find_packages

setup(
    name='wopmetabarcoding',
    version='0.1',
    license='MIT',
    author='Aitor Gonzalez, Thomas Dechatre',
    author_email='aitor.gonzalez@univ-amu.fr',
    url='http://www.aitorgonzalezlab.org',
    long_description="readme.md",
    packages=find_packages(),
    include_package_data=True,
    description="Metabarcoding wrappers and models for WopMars",
    install_requires=['biopython'],
    )

