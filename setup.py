from setuptools import setup
from setuptools import find_packages

install_requires = [
		'biopython>=1.14.2',
		'wopmars>=1.1.29',
		'Jinja2>=2.10',
		'pandas==0.23.3',
]

setup(
    name='wopmetabarcoding',
    version='0.2.0',
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
        'console_scripts': ['create_db_taxonomy=bin.create_db_taxonomy:main']
    },
    )

