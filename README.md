% VTAM - 

# Installation

## Build environment setup with Conda

Create a conda environment with vsearch and blast

~~~
conda create --name vtam python=3.7
~~~

Activate conda

~~~
conda activate vtam
~~~

## Install dependencies and VTAM

List of dependencies

- vsearch=2.7.0
- blast=2.9.0
- wopmars=0.0.8 ( https://github.com/aitgon/wopmars )

These dependencies and VTAM can be installed automatically using within the Conda environment:

~~~
make
~~~

# Documentation (TODO)

Python docstrings use the google style.

To build the documentation, run these commands

~~~
cd doc
sphinx-apidoc -o . ../vtam
make html
~~~

