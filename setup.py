# -*- coding: UTF-8 -*-

__author__ = "Aitor Gonzalez, Thomas Dechatre, Reda Mekdad, Emese Meglecz"
__copyright__ = "Copyright 2018-2020, Aitor Gonzalez, Emese Meglecz"
__email__ = "aitor.gonzalez@univ-amu.fr, emese.meglecz@univ-amu.fr"
__license__ = "MIT"

import configparser
from setuptools import setup
from setuptools import find_packages
import os
import sys

config = configparser.RawConfigParser()
config.read(os.path.join('.', 'setup.cfg'))
version = config['metadata']['version']
author = config['metadata']['author']
email = config['metadata']['email']
license = config['metadata']['license']

if sys.version_info < (3, 6):
    print("At least Python 3.6 is required.\n", file=sys.stderr)
    exit(1)

try:
    from setuptools import setup, find_packages
except ImportError:
    print("Please install setuptools before installing VTAM.",
          file=sys.stderr)
    exit(1)

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as fin:
    long_description = fin.read()

CLASSIFIERS = """\
Development Status :: 4 - Beta
Environment :: Console
Intended Audience :: Science/Research
License :: OSI Approved :: MIT License
Programming Language :: Python
Programming Language :: Python :: 3
Programming Language :: Python :: 3.7
Programming Language :: Python :: 3 :: Only
Topic :: Scientific/Engineering :: Bio-Informatics
Operating System :: POSIX :: Linux
Operating System :: Microsoft :: Windows :: Windows 10
"""


# Create list of package data files
def data_files_to_list(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths


data_file_list = data_files_to_list('vtam/data')

setup(
    name='vtam',
    version=version,
    description="VTAM - Validation and Taxonomic Assignation of Metabarcoding Data",
    author=author,
    author_email=email,
    url="https://vtam.readthedocs.io/en/latest/",
    license=license,
    long_description=long_description,
    classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
    packages=find_packages(),
    package_dir={'vtam': 'vtam'},
    package_data={'vtam': data_file_list},
    install_requires=['jinja2', 'pyyaml', 'sqlalchemy', 'biopython', 'pandas', 'termcolor', 'wopmars'],
    entry_points={
        'console_scripts': ['vtam=vtam:main']
    },
)
