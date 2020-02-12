VTAM - Validation and Taxonomic Assignation of Metabarcoding Data
=================================================================

Installation
------------

See the 'tutorial/tutorial.md'

Documentation
--------------

Python docstrings use the google style.

To build the documentation, you need jinja2, sphinx and sphinx-rtd-theme

Then, run these commands

.. code:: bash

    cd doc
    sphinx-apidoc -o . ../vtam
    make html

