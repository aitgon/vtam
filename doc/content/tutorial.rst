Tutorial
============

Installation
---------------------------------------------------------------------------------

The installation instructions can be found in the readme file.

Data
---------------------------------------------------------------------------------

To setup a new VTAM workflow, follow these steps:

Prepare a FASTQ directory with the FASTQ dir. Here we will use a dataset from our previous publication: `PMID 28776936 <https://pubmed.ncbi.nlm.nih.gov/28776936>`_.

The dataset can be found as Dryad dataset ` doi:10.5061/dryad.f40v5 <https://datadryad.org/stash/dataset/doi:10.5061/dryad.f40v5>`_.

1. Set FASTQ file location
2. Define FASTQ file sample information

Merge the FASTQ files
---------------------------------------------------------------------------------

Set FASTQ file location
^^^^^^^^^^^^^^^^^^^^^^^^

FASTQ files are in this directory

.. code-block:: bash

    $ ls /path/to/my/fastqdir
    MFZR1_S4_L001_R1_001.fastq  MFZR2_S5_L001_R2_001.fastq  ZFZR1_S1_L001_R1_001.fastq  ZFZR2_S2_L001_R2_001.fastq
    MFZR1_S4_L001_R2_001.fastq  MFZR3_S6_L001_R1_001.fastq  ZFZR1_S1_L001_R2_001.fastq  ZFZR3_S3_L001_R1_001.fastq
    MFZR2_S5_L001_R1_001.fastq  MFZR3_S6_L001_R2_001.fastq  ZFZR2_S2_L001_R1_001.fastq  ZFZR3_S3_L001_R2_001.fastq

Define FASTQ file sample information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a TSV (tab-separated file), with a header and 10 columns with all the information per FASTQ file pair.

These columns are needed

- Tag_fwd
- Primer_fwd
- Tag_rev
- Primer_rev
- Marker
- Biosample
- Replicate
- Run
- Fastq_fw
- Fastq_rv


The first two lines of my *fastqinfo.tsv* look like this:

.. code-block:: bash

    Tag_fwd	Primer_fwd	Tag_rev	Primer_rev	Marker	 Biosample	Replicate	Run	Fastq_fwd	Fastq_rev
    cgatcgtcatcacg	TCCACTAATCACAARGATATTGGTAC	cgcgatctgtagag	WACTAATCAATTWCCAAATCCTCC	MFZR	14Mon01	repl2	prerun	MFZR2_S5_L001_R1_001.fastq	MFZR2_S5_L001_R2_001.fastq

Run the VTAM merge command
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to */path/to/my/fastqdir* and *fastqinfo.tsv*, we need

- Output WopMars DB file
- Output TSV file with FASTA file sample information
- Output directory to write the merged FASTA files

.. code-block:: bash

    mkdir -p out/fastadir
    vtam merge --fastqinfo fastqinfo.tsv --fastqdir /path/to/my/fastqdir --fastainfo out/fastainfo.tsv --fastadir out/fastadir --log out/vtam.log -v


Open the *fastainfo.tsv* file and verify its content. A new column should be written with the names of the merged FASTA files.

Verify also the content of the *out/fastadir* with the merged FASTA files.

Demultiplex the reads, filter variants and create the ASV tables
-------------------------------------------------------------------------

There is a single command *filter* to demultiplex the reads, filter variants and create the ASV tables. This command takes quite long but its progress can be seen in the log file.

.. code-block:: bash

    vtam filter --fastainfo out/fastainfo.tsv --fastadir out/fastadir --db out/db.sqlite --outdir out --log out/vtam.log -v

The variants that passed all the filters together with read count in the different biosamples are found in the *out/asvtable.tsv*. The variants that were removed by the different filters can be found in the *out/db.sqlite* database that can be opened with the *sqlitebrowser* program.

Pool Markers
----------------

When variants were amplified with different markers, these variants can be pooled around a variant centroid with the following commands.

An input TSV file must be given with the run and marker combinations that must be pooled. Eg, this is the *pool_run_marker.tsv* file:

.. code-block:: bash

    run_name	marker_name
    prerun	MFZR
    prerun	ZFZR

Then the *pool_markers* subcommand can be used:


.. code-block:: bash

    vtam pool_markers --db ${DB} --runmarker pool_run_marker.tsv --pooledmarkers out/pooled_markers.tsv

Taxon Assignation
-------------------------------------------------------------------------

There is the 'taxassign' subcommand that can assign taxa. 

To assign variants to taxa, we need the COI blast DB and the taxonomy information.

Precomputed versions of these files can be downloaded in the following way:

.. code-block:: bash

    vtam taxonomy -o out/taxonomy.tsv --precomputed
    vtam coi_blast_db --coi_blast_db_dir out/coi_blast_db

The input file is a TSV file, where the last column are the sequence of the variants. Both the *out/asvtable.tsv* and *pool_run_marker.tsv* can be used for the assignation.

The command to carry out the taxon assignation is:

.. code-block:: bash

    vtam taxassign --variants out/pooled_markers.tsv --variant_taxa out/pooled_markers_taxa.tsv --db out/db.sqlite --taxonomy out/taxonomy.tsv --blastdbdir out/coi_blast_db --blastdbname coi_blast_db --log out/vtam.log

Parameter Optimization
------------------------------

To help the user select the parameters, VTAM has an *optimize* subcommand that will compute different values based on positive and negative variants present in the mock, negative and real biosamples. The set of known variants are defined in a TSV file like this: :download:`known_occurrences.tsv <known_occurrences.tsv>`

.. code-block:: bash

    vtam optimize --fastainfo out/fastainfo.tsv --fastadir out/fastadir --known_occurrences known_occurrences.tsv --db out/db.sqlite --outdir out --log out/vtam.log -v


