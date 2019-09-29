from setuptools import setup
from setuptools import find_packages

install_requires = [
    'biopython>=1.14.2',
    'wopmars>y=0.0.5',
    'Jinja2>=2.10',
    'pandas==0.23.3',
    'termcolor>=1.1.0',
]

setup(
    name='wopmetabarcoding',
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
        'console_scripts': ['vtam=bin.vtam:main',
                            'vtam_merge=bin.vtam_merge:main',
                            'vtam_optimize=bin.vtam_optimize:main',
                            'create_db_taxonomy=bin.create_db_taxonomy:main',
                            'create_db_accession2taxid=bin.create_db_accession2taxid:main']
    },
    )

