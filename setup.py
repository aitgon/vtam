# -*- coding: UTF-8 -*-

__author__ = "Aitor Gonzalez, Thomas Dechatre, Reda Mekdad, Emese Meglecz"
__copyright__ = "Copyright 2018-2020, Aitor Gonzalez, Emese Meglecz"
__email__ = "aitor.gonzalez@univ-amu.fr, emese.meglecz@univ-amu.fr"
__license__ = "MIT"

from configparser import RawConfigParser
from setuptools import setup
from setuptools import find_packages
import os
import sys

install_requires = [
    'Jinja2>=2.11',
    'PyYAML>=5.3',
    'SQLAlchemy>=1.3',
    'biopython>=1.76',
    'pandas>=1.0',
    'termcolor>=1.1',
    'wopmars>=0.0.11',
]

if sys.version_info < (3, 7):
    print("At least Python 3.7 is required.\n", file=sys.stderr)
    exit(1)

try:
    from setuptools import setup, find_packages
except ImportError:
    print("Please install setuptools before installing wopmars.",
          file=sys.stderr)
    exit(1)

def get_version():
    """Return package version from setup.cfg."""
    config = RawConfigParser()
    config.read(os.path.join('.', 'setup.cfg'))
    return config.get('metadata', 'version')

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as fin:
    long_description = fin.read()

setup(
    name='vtam',
    version=str(get_version()),
    description="VTAM - Validation and Taxonomic Assignation of Metabarcoding Data",
    long_description=long_description,
    license=__license__,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    author=__author__,
    author_email=__email__,
    url='https://tagc.univ-amu.fr/en/users/gonzalez-aitor, http://net.imbe.fr/~emeglecz/',
    packages=find_packages(),
    package_dir={'vtam': 'vtam'},
    package_data={'vtam': ["data/*.yml", "tests/*.py", "tests/test_files/*", "tests/test_files/optimize_f7/*",
                           "tests/output/*"]},
    data_files = [],
    install_requires=install_requires,
    entry_points={
        'console_scripts': ['vtam=vtam:main']
    },
)
