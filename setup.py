# -*- coding: UTF-8 -*-

__author__ = "Aitor Gonzalez, Thomas Dechatre, Reda Mekdad, Emese Meglecz"
__copyright__ = "Copyright 2018-2020, Aitor Gonzalez, Emese Meglecz"
__email__ = "aitor.gonzalez@univ-amu.fr, emese.meglecz@univ-amu.fr"
__license__ = "MIT"

from configparser import RawConfigParser
from setuptools import setup
from setuptools import find_packages
# from vtam import __version__
import os
import sys
import vtam


def read_setup_cfg_metadata(field):
    """Return package version from setup.cfg."""
    config = RawConfigParser()
    config.read(os.path.join('.', 'setup.cfg'))
    return str(config.get('metadata', field))

install_requires = [
    'Jinja2>=2.11',
    'PyYAML>=5.3',
    'SQLAlchemy>=1.3',
    'biopython>=1.76',
    'pandas>=1.0',
    'termcolor>=1.1',
    'wopmars>=0.1.1',
]

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
    # version=read_setup_cfg_metadata(field='version'),
    version=vtam.__version__,
    description="VTAM - Validation and Taxonomic Assignation of Metabarcoding Data",
    author=read_setup_cfg_metadata(field='author'),
    author_email=read_setup_cfg_metadata(field='email'),
    url='https://tagc.univ-amu.fr/en/users/gonzalez-aitor, http://net.imbe.fr/~emeglecz/',
    license=read_setup_cfg_metadata(field='license'),
    long_description=long_description,
    classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
    packages=find_packages(),
    package_dir={'vtam': 'vtam'},
    package_data={'vtam': data_file_list},
    install_requires=install_requires,
    entry_points={
        'console_scripts': ['vtam=vtam:main']
    },
)
