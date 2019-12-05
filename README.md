# Documentation

Python docstrings use the google style.

To build the documentation, run these commands

~~~
cd doc
sphinx-apidoc -o . ../vtam
make html
~~~

# Install

~~~
conda deactivate
conda create --name vtam python=3.7
conda activate vtam
conda install -c bioconda vsearch=2.7.0
conda install -c bioconda blast=2.9.0

pip install -e .

make
~~~


