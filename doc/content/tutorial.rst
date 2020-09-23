Tutorial
============

<<<<<<< HEAD
Important! Before running any command, do not forget to change directory to vtam-X.Y.Z and activate the conda environment.

.. code-block:: bash

    cd vtam-X.Y.Z
    conda activate vtam

.. note::
    With the exception of BLAST database files and the sqlite database all I/O files of VTAM are text files, that can be opened and edited by a simple text editor (gedit, geany, Notepad++ etc.):

    - TSV:  Text files with tab separated values. Can also be opened by spreadsheets such as LibreOffice, Excel
    - YML: Text files used to provide parameter names and values

Data
---------------------------------------------------------------------------------

In this tutorial, we use a small test dataset based on our previous publication: `PMID 28776936 <https://pubmed.ncbi.nlm.nih.gov/28776936>`_. In this dataset, each sample was amplified by two overlapping markers (mfzr and zfzr), targeting the first 175-181 nucleotides of the COI gene (Fig. 2). We had three PCR replicates for each sample-marker combination (Fig. 1). The samples are tagged <link to glossary>, so the combination of the forward and reverse tags can be used to identify the origin (sample) of each read. set, each sample was amplified by two overlapping markers (mfzr and zfzr), targeting the first 175-181
=======
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
>>>>>>> 2195a5d6bb2972454c238c6fd43a8da854059011

To reduce run time, the test dataset contains only one mock sample, one negative control and two real samples.

.. figure:: img/tuto_fig2.png
   :scale: 50 %
   :alt: Figure 2

<<<<<<< HEAD
Figure 2. Positions on the primer on the COI gene used in the test dataset.

You can download these FASTQ files from here with this command:

.. code-block:: bash

    wget -nc http://pedagogix-tagc.univ-mrs.fr/~gonzalez/vtam/fastq.tar.gz -O fastq.tar.gz
    tar zxvf fastq.tar.gz
    rm fastq.tar.gz

This will create a "FASTQ" directory with 12 FASTQ files:

.. code-block:: bash

    $ ls fastq/
    mfzr_1_fw.fastq
    mfzr_1_rv.fastq
    ...

mfzr_1_fw.fastq: Forward reads of replicates of the MFZR marker (all samples)

merge: Merge FASTQ files
----------------------------

The simplest use of vtam is to analyze one sequencing run (run1) and one marker (MFZR).
=======
    $ ls /path/to/my/fastqdir
    MFZR1_S4_L001_R1_001.fastq  MFZR2_S5_L001_R2_001.fastq  ZFZR1_S1_L001_R1_001.fastq  ZFZR2_S2_L001_R2_001.fastq
    MFZR1_S4_L001_R2_001.fastq  MFZR3_S6_L001_R1_001.fastq  ZFZR1_S1_L001_R2_001.fastq  ZFZR3_S3_L001_R1_001.fastq
    MFZR2_S5_L001_R1_001.fastq  MFZR3_S6_L001_R2_001.fastq  ZFZR2_S2_L001_R1_001.fastq  ZFZR3_S3_L001_R2_001.fastq

Define FASTQ file sample information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
>>>>>>> 2195a5d6bb2972454c238c6fd43a8da854059011

The first step is to merge <link to glossary> the FASTQ files and transform them into fasta files. It can be skipped, if you have single end reads, or your paired sequences have already been merged and transformed into fasta files.

Create a TSV (tab-separated file), with a header and 10 columns with all the information per FASTQ file pair. We will call it "fastqinfo_mfzr.tsv" in this tutorial and you will find it in the VTAM "doc/data" directory [Link to file in "doc directory"]. This TSV file will determine, which file pairs should be merged. These files should be all in the "fastq" directory. This directory can contain other files as well, but they will not be analyzed.

The following columns are required:

- TagFwd
- PrimerFwd
- TagRev
- PrimerRev
- Marker
- Sample
- Replicate
- Run
<<<<<<< HEAD
- FastqFwd
- FastqRev

Tag and primer sequences are in 5' => 3' orientation.

Hereafter are the first lines of the "fastqinfo_mfzr.tsv" file:

.. code-block:: bash

    TagFwd    PrimerFwd    TagRev    PrimerRev    Marker    Sample    Replicate    Run    FastqFwd    FastqRev
    tcgatcacgatgt    TCCACTAATCACAARGATATTGGTAC    tgtcgatctacagc    WACTAATCAATTWCCAAATCCTCC    mfzr    tpos1_run1    1    run1    mfzr_1_fw.fastq    mfzr_1_rv.fastq
    agatcgtactagct    TCCACTAATCACAARGATATTGGTAC    tgtcgatctacagc    WACTAATCAATTWCCAAATCCTCC    mfzr    tnegtag_run1    1    run1    mfzr_1_fw.fastq    mfzr_1_rv.fastq

We propose to work in a project directory called "asper1" (the dataset comes from a project on Zingel asper) and copy user created input files such as "fastqinfo_mfzr.tsv" to the "asper1/user_input" directory.

.. code-block:: bash

    asper1
    `-- user_input
      `-- fastqinfo_mfzr.tsv
    fastq
    |-- mfzr_1_fw.fastq
    |-- mfzr_1_rv.fastq
    |-- ...

Run merge for all file-pairs in the "fastqinfo_mfzr.tsv"

.. code-block:: bash

    vtam merge --fastqinfo asper1/user_input/fastqinfo_mfzr.tsv --fastqdir fastq --fastainfo asper1/run1_mfzr/fastainfo.tsv --fastadir asper1/run1_mfzr/merged -v --log asper1/vtam.log

.. note::
    For info on I/O files see the Reference section<LINK The merge command>

This command adds a "merged" directory and a new "fastainfo_mfzr.tsv" file:

.. code-block:: bash

    asper1
    |-- run1_mfzr
    |  |-- fastainfo.tsv
    |  `-- merged
    |    |-- mfzr_1_fw.fasta
    |    |-- mfzr_2_fw.fasta
    |    `-- mfzr_3_fw.fasta
    |-- user_input
    |  |-- fastqinfo_mfzr.tsv
    |-- vtam.err
    `-- vtam.log
    fastq
    |-- mfzr_1_fw.fastq
    |-- mfzr_1_rv.fastq
    |-- ...

The first lines of the "fastainfo_mfzr.tsv" look like this:

    run    marker    sample    replicate    tagfwd    primerfwd    tagrev    primerrev    mergedfasta
    run1    mfzr    tpos1_run1    1    tcgatcacgatgt    TCCACTAATCACAARGATATTGGTAC    tgtcgatctacagc    WACTAATCAATTWCCAAATCCTCC    mfzr_1_fw.fasta
    run1    mfzr    tnegtag_run1    1    agatcgtactagct    TCCACTAATCACAARGATATTGGTAC    tgtcgatctacagc    WACTAATCAATTWCCAAATCCTCC    mfzr_1_fw.fasta

sortreads: Demultiplex and trim the reads
--------------------------------------------------------

There is a single command "sortreads" to demultiplex <link to glossary> the reads according to tags and to trim off tags and primers.

The sortreads command is designed to deal with a dual indexing, where forward and reverse tag combinations are used to determine the origin of the reads. This is one of the most complex case of demultiplexing, therefore we implemented "sortreads" to help users.

For simpler cases, we suggest using "cutadapt" directly, since it is quite straightforward.

.. code-block:: bash

    vtam sortreads --fastainfo asper1/run1_mfzr/fastainfo.tsv --fastadir asper1/run1_mfzr/merged --sorteddir asper1/run1_mfzr/sorted -v --log asper1/vtam.log

.. note::
    For info on I/O files see the Reference section<LINK The sortreads command>

The FASTA files with the sorted reads are written to the "asper1/sorted" directory:

.. code-block:: bash

    asper1
    |-- run1_mfzr
    |  |-- fastainfo.tsv
    |  |-- ...
    |  `-- sorted
    |    |-- mfzr_1_fw_000.fasta
    |    |-- mfzr_1_fw_001.fasta
    |    |-- ...
    |    `-- sortedinfo.tsv
    |-- ...
    ...

In addition, the TSV file "asper1/run1_mfzr/sorted/sortedinfo.tsv" lists the information, i.e. run, marker, sample and replicate about each sorted FASTA file. The "sortedinfo.tsv" file looks like this:

.. code-block:: bash

    run    marker    sample    replicate    sortedfasta
    run1    MFZR    tpos1_run1    1    mfzr_1_fw_000.fasta
    run1    MFZR    tnegtag_run1    1    mfzr_1_fw_001.fasta

filter: Filter variants and create the ASV table
---------------------------------------------------
=======
- Fastq_fw
- Fastq_rv


The first two lines of my *fastqinfo.tsv* look like this:

.. code-block:: bash

    Tag_fwd	Primer_fwd	Tag_rev	Primer_rev	Marker	 Sample	Replicate	Run	Fastq_fwd	Fastq_rev
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
>>>>>>> 2195a5d6bb2972454c238c6fd43a8da854059011

The "filter" command is typically first run with default parameters. From the output, users should identify clearly unwanted (‘delete’) and clearly necessary (‘keep’) occurrences (see Manual section for details). These false positive and false negative occurrences will be used as input to the "optimize" command. The "optimize" command will then suggest an optimal parameter combination tailored to your dataset. Then "filter" command should be run again with the optimized parameters.

<<<<<<< HEAD
Let's run first the "filter" command with default parameters.
=======
Open the *fastainfo.tsv* file and verify its content. A new column should be written with the names of the merged FASTA files.
>>>>>>> 2195a5d6bb2972454c238c6fd43a8da854059011

Verify also the content of the *out/fastadir* with the merged FASTA files.

<<<<<<< HEAD
    vtam filter --db asper1/db.sqlite --sortedinfo asper1/run1_mfzr/sorted/sortedinfo.tsv --sorteddir asper1/run1_mfzr/sorted --asvtable asper1/run1_mfzr/asvtable_default.tsv -v --log asper1/vtam.log

.. note::
    For info on I/O files see the Reference section<LINK The filter command>

This command creates two new files "db.sqlite" and "asvtable_mfzr_default.tsv":

.. code-block:: bash

    asper1
    |-- db.sqlite
    |-- run1_mfzr
    |  |-- asvtable_default.tsv
    |-- ...
    ...

The database "asper1/db.sqlite" contains one table by filter, and in each table occurrences are marked as deleted (filter_delete = 1) or retained  (filter_delete = 0). This database can be opened with a sqlite browser program (For example, https://sqlitebrowser.org / or https://sqlitestudio.pl/index.rvt).

.. figure:: img/tuto_fig3.png
   :scale: 50 %
   :alt: Figure 3

The "asper1/run1_mfzr/asvtable_default.tsv" contains information about the variants that passed all the filters such as the run, maker, read count over all replicates of a sample and the sequence. Hereafter are the first lines of the asvtable_default.tsv
=======
Demultiplex the reads, filter variants and create the ASV tables
-------------------------------------------------------------------------

There is a single command *filter* to demultiplex the reads, filter variants and create the ASV tables. This command takes quite long but its progress can be seen in the log file.

.. code-block:: bash

    vtam filter --fastainfo out/fastainfo.tsv --fastadir out/fastadir --db out/db.sqlite --outdir out --log out/vtam.log -v

The variants that passed all the filters together with read count in the different samples are found in the *out/asvtable.tsv*. The variants that were removed by the different filters can be found in the *out/db.sqlite* database that can be opened with the *sqlitebrowser* program.
>>>>>>> 2195a5d6bb2972454c238c6fd43a8da854059011

Pool Markers
----------------

<<<<<<< HEAD
    run    marker    variant    sequence_length    read_count    tpos1_run1    tnegtag_run1    14ben01    14ben02    clusterid    clustersize    chimera_borderline    sequence
    run1    MFZR    25    181    478    478    0    0    0    25    1    False    ACTATACCTTATCTTCGCAGTATTCTCAGGAATGCTAGGAACTGCTTTTAGTGTTCTTATTCGAATGGAACTAACATCTCCAGGTGTACAATACCTACAGGGAAACCACCAACTTTACAATGTAATCATTACAGCTCACGCATTCCTAATGATCTTTTTCATGGTTATGCCAGGACTTGTT
    run1    MFZR    51    181    165    0    0    0    165    51    1    False    ACTATATTTAATTTTTGCTGCAATTTCTGGTGTAGCAGGAACTACGCTTTCATTGTTTATTAGAGCTACATTAGCGACACCAAATTCTGGTGTTTTAGATTATAATTACCATTTGTATAATGTTATAGTTACGGGTCATGCTTTTTTGATGATCTTTTTTTTAGTAATGCCTGCTTTATTG


.. note::
    Filter can be run with the --known_occurrences argument that will add an additional column for each mock sample flagging expected variants. This helps in creating the known_occurrences.tsv input file for the optimization step <LINK To  Make the ASV table>

taxassign: Assign variants of ASV table to taxa
--------------------------------------------------

The "taxassign" command assigns ASV sequences in the last column of a TSV file such as the "asvtable_default.tsv" file to taxa.

The "taxassign" command needs a BLAST database (containing reference sequences of known taxonomic origin) <LINK to creation of taxonomic database> and the taxonomy information file <LINK to The command taxonomy and the taxonomic lineage input>.

A precomputed taxonomy file in TSV format and the BLAST database with COI sequences can be downloaded with these commands:

.. code-block:: bash

    vtam taxonomy -output vtam_db/taxonomy.tsv --precomputed
    vtam coi_blast_db --blastdbdir vtam_db/coi_blast_db

These commands result in these new files:

.. code-block:: bash

    ...
    vtam_db
    |-- coi_blast_db
    |  |-- coi_blast_db_20200420.nhr
    |  |-- coi_blast_db_20200420.nin
    |  |-- coi_blast_db_20200420.nog
    |  |-- coi_blast_db_20200420.nsd
    |  |-- coi_blast_db_20200420.nsi
    |  └-- coi_blast_db_20200420.nsq
    `-- taxonomy.tsv

.. note::
    Alternatively, you can use your own custom database or the NCBI nucleotide database <LINK to the tax database creation>

Then, we can carry out the taxonomic assignation of variants in the "asvtable_default.tsv" with the following command:

.. code-block:: bash
=======
When variants were amplified with different markers, these variants can be pooled around a variant centroid with the following commands.

An input TSV file must be given with the run and marker combinations that must be pooled. Eg, this is the *pool_run_marker.tsv* file:

.. code-block:: bash

    run_name	marker_name
    prerun	MFZR
    prerun	ZFZR

Then the *pool_markers* subcommand can be used:
>>>>>>> 2195a5d6bb2972454c238c6fd43a8da854059011

    vtam taxassign --db asper1/db.sqlite --asvtable asper1/run1_mfzr/asvtable_default.tsv --output asper1/run1_mfzr/asvtable_default_taxa.tsv --taxonomy vtam_db/taxonomy.tsv --blastdbdir vtam_db/coi_blast_db --blastdbname coi_blast_db_20200420 -v --log asper1/vtam.log

<<<<<<< HEAD
.. note::
    For info on I/O files see the Reference section<LINK The taxassign command>

This results in an additional file:

.. code-block:: bash

    asper1/
    |-- run1_mfzr
    |  |-- asvtable_default.tsv
    |  |-- asvtable_default_taxa.tsv

optimize: Compute optimal filter parameters based on mock and negative samples
---------------------------------------------------------------------------------------

The "optimize" command helps users choose optimal parameters for filtering that are specifically adjusted to the dataset. This optimization is based on mock samples and negative controls.

Users should prepare a TSV file ("known_occurences_mfzr.tsv") with occurrences to be kept in the results (typically expected variants of the mock samples) and occurrences to be clearly deleted (typically all occurrences in negative controls, and unexpected occurrences in the mock samples).<link to reference section of the Manual>

The example TSV file for the known occurrences of the MFZR marker can be found in the "doc/data" directory <link>.

The first lines of this file look like this:

.. code-block:: bash

    Marker    Run    Sample    Mock    Variant    Action    Sequence
    MFZR    run1    tpos1_run1    1        keep    ACTATATTTTATTTTTGGGGCTTGATCCGGAATGCTGGGCACCTCTCTAAGCCTTCTAATTCGTGCCGAGCTGGGGCACCCGGGTTCTTTAATTGGCGACGATCAAATTTACAATGTAATCGTCACAGCCCATGCTTTTATTATGATTTTTTTCATGGTTATGCCTATTATAATC
    MFZR    run1    tpos1_run1    1        keep    ACTTTATTTTATTTTTGGTGCTTGATCAGGAATAGTAGGAACTTCTTTAAGAATTCTAATTCGAGCTGAATTAGGTCATGCCGGTTCATTAATTGGAGATGATCAAATTTATAATGTAATTGTAACTGCTCATGCTTTTGTAATAATTTTCTTTATAGTTATACCTATTTTAATT

    ...
    MFZR    run1    tpos1_run1    1        delete    TTTATATTTCATTTTTGGTGCATGATCAGGTATGGTGGGTACTTCCCTTAGTTTATTAATTCGAGCAGAACTTGGTAATCCTGGTTCTTTGATTGGCGATGATCAGATTTATAACGTTATTGTCACTGCCCATGCTTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT
    MFZR    run1    tnegtag_run1    0        delete    TTTATATTTTATTTTTGGAGCCTGAGCTGGAATAGTAGGTACTTCCCTTAGTATACTTATTCGAGCCGAATTAGGACACCCAGGCTCTCTAATTGGAGACGACCAAATTTATAATGTAATTGTTACTGCTCATGCTTTTGTAATAATTTTTTTTATAGTTATGCCAATTATAATT

.. note::

    It is possible to add extra columns with your notes (for example taxon names) to this file after the “Sequence” column. They will be ignored by VTAM.

The "optimize" command is run like this:

.. code-block:: bash

    vtam optimize --db asper1/db.sqlite --sortedinfo asper1/run1_mfzr/sorted/sortedinfo.tsv --sorteddir asper1/run1_mfzr/sorted --known_occurrences asper1/user_input/known_occurrences_mfzr.tsv --outdir asper1/run1_mfzr -v --log asper1/vtam.log

.. note::

    For info on I/O files see the Reference section<LINK The optimize command>

This command creates four new files:

.. code-block:: bash

    asper1/
    |-- db.sqlite
    |-- run1_mfzr
    |  |-- ...
    |  |-- optimize_lfn_sample_replicate.tsv
    |  |-- optimize_lfn_read_count_and_lfn_variant.tsv
    |  |-- optimize_lfn_variant_specific.tsv
    |  |-- optimize_pcr_error.tsv

.. note::

    Running vtam optimize will run three underlying scripts:

    - OptimizePCRerror, to optimize “pcr_error_var_prop”
    - OptimizeLFNsampleReplicate, to optimize “lfn_sample_replicate_cutoff”
    - OptimizeLFNreadCountAndLFNvariant, to optimize “lfn_read_count_cutoff” and “lfn_variant_cutoff”.

While OptimizePCRerror and OptimizeLFNsampleReplicate do not depend on the other two parameters to be optimized, OptimizeLFNreadCountAndLFNvariant does. For a finer tuning, it is possible to run the three subscripts one by one, and use the optimized values of “pcr_error_var_prop” and “lfn_sample_replicate_cutoff” instead of their default values, when running OptimizeLFNreadCountAndLFNvariant. This procedure can propose less stringent values for “lfn_read_count_cutoff” and “lfn_variant_cutoff”, but still eliminate as many as possible unexpected occurrences, and keep all expected ones.

To run just one subscript, the --until flag can be added to the vtam optimize command

- --until OptimizePCRerror
- --unlit OptimizeLFNsampleReplicate
- --until OptimizeLFNreadCountAndLFNvariant

e.g.

.. code-block:: bash

    vtam optimize --db asper1/db.sqlite --sortedinfo asper1/run1_mfzr/sorted/sortedinfo.tsv --sorteddir asper1/run1_mfzr/sorted --known_occurrences asper1/user_input/known_occurrences_mfzr.tsv --outdir asper1/run1_mfzr -v --log asper1/vtam.log --until OptimizePCRerror

    vtam optimize --db asper1/db.sqlite --sortedinfo asper1/run1_mfzr/sorted/sortedinfo.tsv --sorteddir asper1/run1_mfzr/sorted --known_occurrences asper1/user_input/known_occurrences_mfzr.tsv --outdir asper1/run1_mfzr -v --log asper1/vtam.log --until OptimizeLFNsampleReplicate

Create a params_optimize_mfzr.yml file that will contain the optimal values suggested for “lfn_sample_replicate_cutoff” and “pcr_error_var_prop”
lfn_sample_replicate_cutoff: 0.003
pcr_error_var_prop: 0.1

Run OptimizeLFNreadCountAndLFNvariant with the optimized parameters for the above two parameters.

.. code-block:: bash

    vtam optimize --db asper1/db.sqlite --sortedinfo asper1/run1_mfzr/sorted/sortedinfo.tsv --sorteddir asper1/run1_mfzr/sorted --known_occurrences asper1/user_input/known_occurrences_mfzr.tsv --outdir asper1/run1_mfzr -v --log asper1/vtam.log --until OptimizeLFNreadCountAndLFNvariant --params asper1/user_input/params_optimize_mfzr.yml

This step will suggest the following parameter values
lfn_variant_cutoff: 0.001
lfn_read_count_cutoff: 20
For simplicity, we continue the tutorial with parameters optimized previously, with running all 3 optimize steps in one command.

filter: Create an ASV table with optimal parameters and assign variants to taxa
---------------------------------------------------------------------------------

Once the optimal filtering parameters are chosen, rerun the "filter" command using the existing "db.sqlite" database that already has all the variant counts.

Make a "params_mfzr.yml" file that contains the parameter names and values that differ from the default settings.

The "params_mfzr.yml" can be found in the "doc/data" directory and it looks like this:

.. code-block:: bash

    lfn_variant_cutoff: 0.001
    lfn_sample_replicate_cutoff: 0.003
    lfn_read_count_cutoff: 70
    pcr_error_var_prop: 0.1

Run filter with optimized parameters:

.. code-block:: bash

    vtam filter --db asper1/db.sqlite --sortedinfo asper1/run1_mfzr/sorted/sortedinfo.tsv --sorteddir asper1/run1_mfzr/sorted --params asper1/user_input/params_mfzr.yml --asvtable asper1/run1_mfzr/asvtable_optimized.tsv -v --log asper1/vtam.log

Running again "taxassign" will complete the "asvtable_optimized.tsv" with the taxonomic information. It will be very quick since most variants in the table have already gone through the taxonomic assignment, and these assignations are extracted from the “db.sqlite”.

.. code-block:: bash

    vtam taxassign --db asper1/db.sqlite --asvtable asper1/run1_mfzr/asvtable_optimized.tsv --output asper1/run1_mfzr/asvtable_optimized_taxa.tsv --taxonomy vtam_db/taxonomy.tsv --blastdbdir vtam_db/coi_blast_db --blastdbname coi_blast_db_20200420 -v --log asper1/vtam.log

We finished our first analysis with VTAM! The resulting directory structure looks like this:

.. code-block:: bash

    asper1/
    |-- db.sqlite
    |-- run1_mfzr
    |  |-- asvtable_default.tsv
    |  |-- asvtable_default_taxa.tsv
    |  |-- asvtable_optimized.tsv
    |  |-- asvtable_optimized_taxa.tsv
    |  |-- fastainfo.tsv
    |  |-- merged
    |  |  |-- mfzr_1_fw.fasta
    |  |  |-- ...
    |  |-- optimize_lfn_sample_replicate.tsv
    |  |-- optimize_lfn_read_count_and_lfn_variant.tsv
    |  |-- optimize_lfn_variant_specific.tsv
    |  |-- optimize_pcr_error.tsv
    |  `-- sorted
    |    |-- mfzr_1_fw_000.fasta
    |    |-- ...
    |    `-- sortedinfo.tsv

Add new run-marker data to the existing database
-----------------------------------------------------

The same samples can be amplified by different but strongly overlapping markers. In this case, it makes sense to pool all the data into the same database, and produce just one ASV table, with information of both markers. This is the case in our test dataset.

It is also frequent to have different sequencing runs (with one or several markers) that are part of the same study. Feeding them to the same database assures coherence in variant IDs, and gives the possibility to easily produce one ASV table with all the runs and avoids re-running the taxassign on variants that have already been assigned to a taxon.

**We assume that you have gone through the basic pipeline in the previous section** [LINK].
Let's see an example on how to complete the previous analyses with the dataset obtained for the same samples but for another marker (ZFZR). The principle is the same if you want to complete the analyses with data from a different sequencing run.
First we need to prepare these user inputs:
The directory with the FASTQ files: "fastqinfo_zfzr.tsv"

This is the "merge" command for the new run-marker:

.. code-block:: bash

    vtam merge --fastqinfo asper1/user_input/fastqinfo_zfzr.tsv --fastqdir fastq --fastainfo asper1/run1_zfzr/fastainfo.tsv --fastadir asper1/run1_zfzr/merged -v --log asper1/vtam.log

This is the "sortreads" command for the new marker ZFZR:

.. code-block:: bash

    vtam sortreads --fastainfo asper1/run1_zfzr/fastainfo.tsv --fastadir asper1/run1_zfzr/merged --sorteddir asper1/run1_zfzr/sorted -v --log asper1/vtam.log

The "filter" command for the new marker ZFZR is the same as in the basic pipeline, but we will complete the previous database "asper1/db.sqlite" with the new variants.

.. code-block:: bash

    vtam filter --db asper1/db.sqlite --sortedinfo asper1/run1_zfzr/sorted/sortedinfo.tsv --sorteddir asper1/run1_zfzr/sorted --asvtable asper1/run1_zfzr/asvtable_default.tsv -v --log asper1/vtam.log

Next we run the "taxassign" command for the new ASV table "asper1/asvtable_zfzr_default.tsv":

.. code-block:: bash

    vtam taxassign --db asper1/db.sqlite --asvtable asper1/run1_zfzr/asvtable_default.tsv --output asper1/run1_zfzr/asvtable_default_taxa.tsv --taxonomy vtam_db/taxonomy.tsv --blastdbdir vtam_db/coi_blast_db --blastdbname coi_blast_db_20200420 -v --log asper1/vtam.log
=======
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

    vtam taxassign --asvtable out/pooled_markers.tsv --variant_taxa out/pooled_markers_taxa.tsv --db out/db.sqlite --taxonomy out/taxonomy.tsv --blastdbdir out/coi_blast_db --blastdbname coi_blast_db --log out/vtam.log

Parameter Optimization
------------------------------

To help the user select the parameters, VTAM has an *optimize* subcommand that will compute different values based on positive and negative variants present in the mock, negative and real samples. The set of known variants are defined in a TSV file like this: :download:`known_occurrences.tsv <known_occurrences.tsv>`

.. code-block:: bash

    vtam optimize --fastainfo out/fastainfo.tsv --fastadir out/fastadir --known_occurrences known_occurrences.tsv --db out/db.sqlite --outdir out --log out/vtam.log -v
>>>>>>> 2195a5d6bb2972454c238c6fd43a8da854059011


