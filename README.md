# Documentation

Python docstrings use the google style.

To build the documentation, run these commands

~~~
cd doc
sphinx-apidoc -o . ../vtam
make html
~~~

# Installation

Create a conda environment with vsearch and blast

~~~
conda deactivate
conda create --name vtam python=3.7
~~~

Activate it

~~~
conda activate vtam
conda install -c bioconda vsearch=2.7.0
conda install -c bioconda blast=2.9.0
~~~

Download and install wopmars

~~~
wget https://github.com/aitgon/wopmars/archive/0.0.8.tar.gz
taz zxvf 0.0.8.tar.gz
pip install wopmars-0.0.8/.
~~~

~~~
pip install .
~~~


