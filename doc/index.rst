.. metabarcoding documentation master file, created by
   sphinx-quickstart on Tue Oct 31 20:17:07 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=============================================================================================
VTAM - Validation and Taxonomic Assignation of Metabarcoding Data
=============================================================================================

.. image:: https://img.shields.io/pypi/v/vtam.svg
    :target: https://pypi.python.org/pypi/vtam

.. image:: https://img.shields.io/pypi/pyversions/vtam.svg
    :target: https://www.python.org

.. image:: https://readthedocs.org/projects/vtam/badge/?version=latest
    :target: http://vtam.readthedocs.io/en/latest/?badge=latest

.. image:: https://github.com/aitgon/vtam/workflows/CI/badge.svg
    :target: https://github.com/aitgon/vtam/actions?query=branch%3Amaster+workflow%3ACI

.. image:: https://travis-ci.org/aitgon/vtam.svg?branch=master
    :target: https://travis-ci.org/aitgon/vtam

.. image:: https://codecov.io/gh/aitgon/vtam/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/aitgon/vtam

VTAM is a metabarcoding package with various commands to process high throughput sequencing (HTS) data of amplicons of one or several metabarcoding markers in FASTQ format and produce a table of amplicon sequence variants (ASVs) assigned to taxonomic groups.
If you use VTAM in scientific works, **please cite the following article**:

**González, A., Dubut, V., Corse, E., Mekdad, R., Dechartre, T. and  Meglécz, E.**. (2020)
`VTAM: A robust pipeline for processing metabarcoding data using internal controls`.
**BioRxiv** . `2020.11.06.371187 <https://www.biorxiv.org/content/10.1101/2020.11.06.371187v3>`_

Commands for a quick installation:

.. code-block:: bash

    conda create --name vtam python=3.7 -y
    python3 -m pip install --upgrade cutadapt
    conda install -c bioconda blast
    conda install -c bioconda vsearch
    python3 -m pip install --upgrade vtam

Commands for a quick working example:

.. code-block:: bash

    vtam example
    cd example
    snakemake --printshellcmds --resources db=1 --snakefile snakefile.yml --cores 4 --configfile asper1/user_input/snakeconfig_mfzr.yml --until asvtable_taxa

The table of amplicon sequence variants (ASV) is here:

.. code-block:: bash

    (vtam) user@host:~/vtam/example$ head -n4 asper1/run1_mfzr/asvtable_default_taxa.tsv
    run	marker	variant	sequence_length	read_count	tpos1_run1	tnegtag_run1	14ben01	14ben02	clusterid	clustersize	chimera_borderlineltg_tax_id	ltg_tax_name	ltg_rank	identity	blast_db	phylum	class	order	family	genus	species	sequence
    run1	MFZR	25	181	478	478	0	0	0	25	1	False	131567	cellular organisms	no rank	80	coi_blast_db_20200420							ACTATACCTTATCTTCGCAGTATTCTCAGGAATGCTAGGAACTGCTTTTAGTGTTCTTATTCGAATGGAACTAACATCTCCAGGTGTACAATACCTACAGGGAAACCACCAACTTTACAATGTAATCATTACAGCTCACGCATTCCTAATGATCTTTTTCATGGTTATGCCAGGACTTGTT
    run1	MFZR	51	181	165	0	0	0	165	51	1	False					coi_blast_db_20200420		ACTATATTTAATTTTTGCTGCAATTTCTGGTGTAGCAGGAACTACGCTTTCATTGTTTATTAGAGCTACATTAGCGACACCAAATTCTGGTGTTTTAGATTATAATTACCATTTGTATAATGTTATAGTTACGGGTCATGCTTTTTTGATGATCTTTTTTTTAGTAATGCCTGCTTTATTG
    run1	MFZR	88	175	640	640	0	0	0	88	1	False	1592914	Caenis pusilla	species	100	coi_blast_db_20200420	Arthropoda	Insecta	Ephemeroptera	Caenidae	Caenis	Caenis pusilla	ACTATATTTTATTTTTGGGGCTTGATCCGGAATGCTGGGCACCTCTCTAAGCCTTCTAATTCGTGCCGAGCTGGGGCACCCGGGTTCTTTAATTGGCGACGATCAAATTTACAATGTAATCGTCACAGCCCATGCTTTTATTATGATTTTTTTCATGGTTATGCCTATTATAATC

The database of intermediate data is here:

.. code-block:: bash

    (vtam) user@host:~/vtam/example$ sqlite3 asper1/db.sqlite '.tables'
   FilterChimera                    Sample
   FilterChimeraBorderline          SampleInformation
   FilterCodonStop                  SortedReadFile
   FilterIndel                      TaxAssign
   FilterLFN                        Variant
   FilterLFNreference               VariantReadCount
   FilterMinReplicateNumber         wom_Execution
   FilterMinReplicateNumber2        wom_FileInputOutputInformation
   FilterMinReplicateNumber3        wom_Option
   FilterPCRerror                   wom_TableInputOutputInformation
   FilterRenkonen                   wom_TableModificationTime
   Marker                           wom_ToolWrapper
   ReadCountAverageOverReplicates   wom_TypeInputOrOutput
   Run

Table of Contents
=================

.. toctree::
    :maxdepth: 2

    content/overview
    content/installation
    content/tutorial
    content/reference
    content/io
    content/glossary
    content/reflist
    content/cite
    content/changelog
    content/contributing

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
