from setuptools import setup
from setuptools import find_packages

install_requires = [
    'biopython>=1.14.2',
    'wopmars>y=0.0.5',
    'Jinja2>=2.10',
    'pandas==0.23.3',
    'SQLAlchemy>=1.3.0',
    'termcolor>=1.1.0',
]

setup(
    name='vtam',
    version='0.0.1',
    license='MIT',
    author='Thomas Dechatre, Aitor Gonzalez',
    author_email='aitor.gonzalez@univ-amu.fr',
    url='http://www.aitorgonzalezlab.org',
    long_description="readme.md",
    packages=find_packages(),
    include_package_data=True,
    description="Metabarcoding wrappers and models for WopMars",
    install_requires=install_requires,
    entry_points={
        'console_scripts': ['vtam=vtam:main',
                            'download_db_blast_coi=vtam.utils.DBblastCOI:DBblastCOI.main',
                            'create_db_taxonomy=vtam.utils.DBtaxonomy:DBtaxonomy.main']
    },
    )

