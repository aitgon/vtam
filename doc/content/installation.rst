Installation
=================================================

Dependencies installation
-------------------------------------------------

Miniconda:

Installation:

.. code-block:: bash

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    conda create --name myenv --python=3.6

Use virtual environment:

.. code-block:: bash

    source activate myenv

    source deactivate

Vscearch:

.. code-block:: bash

    wget https://github.com/torognes/vsearch/archive/v2.8.0.tar.gz
    tar xzf v2.8.0.tar.gz
    cd vsearch-2.8.0
    ./autogen.sh
    ./configure
    make
    make install  # as root or sudo make install

Wopmars:

.. code-block:: bash

    pip install wopmars

Package utilisation
-------------------------------------------------

Step 1: WopFile personalisation

Fill path for input files in the WopFile to launch the
pipeline.

Step 2: Wopmars command

Example:

.. code-block:: bash

    wopmars -w Wopfile.yml -D "sqlite:///data/db.sqlite" -v -p -F

-w: Path to the Wopfile

-D: Path to the database

-v: Set verbosity level

-p: Write logs in standard output

-F: Force the execution of the workflow, without checking for previous executions

In case of problem or for more options write: wopmars -h







