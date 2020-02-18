from setuptools import setup
from setuptools import find_packages

install_requires = [
    'Jinja2>=2.10',
    'PyYAML>=5.1.2',
    'SQLAlchemy>=1.3.0',
    'biopython>=1.14.2',
    'pandas==0.23.3',
    'termcolor>=1.1.0',
    'wopmars>y=0.0.5',
]

setup(
    name='vtam',
    version='0.0.1',
    license='MIT',
    author='Aitor Gonzalez, Thomas Dechatre, Reda Mekdad, Emese Meglecz',
    author_email='aitor.gonzalez@univ-amu.fr',
    url='https://tagc.univ-amu.fr/en/users/gonzalez-aitor',
    long_description="README.rst",
    packages=find_packages(),
    package_dir={'vtam': 'vtam'},
    package_data={'vtam': ["tests/*.py", "tests/test_files/*", "tests/test_files/optimize_f7/*", "tests/output/*"]},  # Add tests
    description="Metabarcoding wrappers and models for WopMars",
    install_requires=install_requires,
    entry_points={
        'console_scripts': ['vtam=vtam:main',]
    },
    )

