=============================================================================================
VTAM - Validation and Taxonomic Assignation of Metabarcoding Data
=============================================================================================

.. image:: https://img.shields.io/pypi/v/vtam.svg?color=blue
    :target: https://pypi.python.org/pypi/vtam

.. image:: https://img.shields.io/pypi/pyversions/vtam.svg
    :target: https://www.python.org

.. image:: https://static.pepy.tech/personalized-badge/vtam?period=month&units=international_system&left_color=gray&right_color=blue&left_text=Downloads
    :target: https://pepy.tech/project/vtam

.. image:: https://codecov.io/gh/aitgon/vtam/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/aitgon/vtam

.. image:: https://readthedocs.org/projects/vtam/badge/?version=latest
    :target: http://vtam.readthedocs.io/en/latest/?badge=latest

.. image:: https://app.travis-ci.com/aitgon/vtam.svg?branch=master
    :target: https://app.travis-ci.com/github/aitgon/vtam

.. image:: https://github.com/aitgon/vtam/workflows/CI/badge.svg
    :target: https://github.com/aitgon/vtam/actions?query=branch%3Amaster+workflow%3ACI

VTAM is a metabarcoding package with various commands to process high throughput sequencing (HTS) data of amplicons of one or several metabarcoding markers in FASTQ format and produce a table of amplicon sequence variants (ASVs) assigned to taxonomic groups.
If you use VTAM in scientific works, **please cite the following article**:

**Aitor González, Vincent Dubut, Emmanuel Corse, Reda Mekdad, Thomas Dechatre, Ulysse Castet, Raphaël Hebert, Emese Meglécz**.
VTAM: A robust pipeline for validating metabarcoding data using controls. Computational and Structural Biotechnology Journal, 2023, 21, pp.1151 - 1156. `10.1016/j.csbj.2023.01.034 <https://dx.doi.org/10.1016/j.csbj.2023.01.034>`_.

For a quick use, the simplest method is to use the `Singularity <https://sylabs.io/singularity>`_ image
that can be downloaded here:
`https://github.com/aitgon/vtam/releases/download/0.2.0/vtam.sif <https://github.com/aitgon/vtam/releases/download/0.2.0/vtam.sif>`_. Singularity can be install as described `Singularity <https://sylabs.io/singularity>`_.

The VTAM container gives access to all commands, e.g.

.. code-block:: bash

    singularity run --app miniconda vtam.sif python --version
    singularity run --app vtam vtam.sif --help
    singularity run --app vtam vtam.sif merge --help

The `VTAM documentation <http://vtam.readthedocs.org/>`_ is hosted at ReadTheDocs.

VTAM is maintained by Aitor González (aitor dot gonzalez at univ-amu dot fr) and Emese Meglécz (emese dot meglecz at univ-amu dot fr).
