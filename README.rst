VTAM - Validation and Taxonomic Assignation of Metabarcoding Data
=================================================================

Installation
------------

Build environment setup with Conda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a conda environment with vsearch and blast

.. code:: bash

    conda create --name vtam python=3.7 -y


Activate conda

.. code:: bash

    conda activate vtam

Fetch and pull latest commit of VTAM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: bash

    git fetch origin master
    git pull origin master

Install dependencies and VTAM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

List of dependencies

* vsearch=2.7.0
* blast=2.9.0
* wopmars=0.0.8 ( https://github.com/aitgon/wopmars )

These dependencies and VTAM can be installed automatically using within the Conda environment:

.. code:: bash

    make

Documentation
--------------

Python docstrings use the google style.

To build the documentation, you need jinja2, sphinx and sphinx-rtd-theme

Then, run these commands

.. code:: bash

    cd doc
    sphinx-apidoc -o . ../vtam
    make html

