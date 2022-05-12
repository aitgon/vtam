# -*- coding: UTF-8 -*-

__author__ = "Aitor González, Raphael Hebert, Thomas Dechatre, Reda Mekdad, Ulysse Castet, Emese Meglecz"
__copyright__ = "Copyright 2018-2022, Aitor González, Emese Meglecz"
__email__ = "aitor.gonzalez@univ-amu.fr, emese.meglecz@univ-amu.fr"
__license__ = "MIT"

import codecs
import configparser
from setuptools import setup
from setuptools import find_packages
import os
import sys

config = configparser.RawConfigParser()
config.read(os.path.join('.', 'setup.cfg'))
author = config['metadata']['author']
email = config['metadata']['email']
license = config['metadata']['license']

cutoff_version_info = (3, 7)
if sys.version_info < cutoff_version_info:
    print("At least Python {}.{} is required.\n".format(cutoff_version_info[0], cutoff_version_info[1]), file=sys.stderr)
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

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")

CLASSIFIERS = """\
Development Status :: 4 - Beta
Environment :: Console
Intended Audience :: Science/Research
License :: OSI Approved :: MIT License
Programming Language :: Python
Programming Language :: Python :: 3
Programming Language :: Python :: 3.7
Programming Language :: Python :: 3.8
Programming Language :: Python :: 3.9
Programming Language :: Python :: 3.10
Programming Language :: Python :: 3 :: Only
Topic :: Scientific/Engineering :: Bio-Informatics
Topic :: Software Development
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
# data_test_list = data_files_to_list('vtam/tests')

setup(
    name='vtam',
    version=get_version("vtam/__init__.py"),
    description="VTAM - Validation and Taxonomic Assignation of Metabarcoding Data "
                "is a metabarcoding pipeline. The analysis starts from high throughput "
                "sequencing (HTS) data of amplicons of one or several metabarcoding "
                "markers and produces an amplicon sequence variant (ASV) "
                "table of validated variants assigned to taxonomic groups.",
    author=author,
    author_email=email,
    url="https://vtam.readthedocs.io",
    license=license,
    long_description=long_description,
    classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
    packages=find_packages(),
    package_dir={'vtam': 'vtam'},
    package_data={'vtam': data_file_list},
    install_requires=['biopython', 'cutadapt', 'jinja2', 'pandas', 'progressbar', 'pyyaml', 'sqlalchemy', 'snakemake', 'termcolor', 'tqdm', 'wopmars'],
    entry_points={
        'console_scripts': ['vtam=vtam:main']
    },
)
