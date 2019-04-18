Tutorial
============

To setup a new VTAM workflow, follow these steps:

1. Set FASTQ file location
2. Define FASTQ file sample information


Set FASTQ file location
---------------------------------------------------------------------------------

FASTQ files are in this directory

.. code-block:: bash

    $ ls /path/to/my/fastqdir
    MFZR1_S4_L001_R1_001.fastq  MFZR2_S5_L001_R2_001.fastq  ZFZR1_S1_L001_R1_001.fastq  ZFZR2_S2_L001_R2_001.fastq
    MFZR1_S4_L001_R2_001.fastq  MFZR3_S6_L001_R1_001.fastq  ZFZR1_S1_L001_R2_001.fastq  ZFZR3_S3_L001_R1_001.fastq
    MFZR2_S5_L001_R1_001.fastq  MFZR3_S6_L001_R2_001.fastq  ZFZR2_S2_L001_R1_001.fastq  ZFZR3_S3_L001_R2_001.fastq


Define FASTQ file sample information
-------------------------------------------------------------------------

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
- Fastq_fwd
- Fastq_rev

The first two lines of my *fastqinfo.tsv* look like this:

.. code-block:: bash

    Tag_fwd	Primer_fwd	Tag_rev	Primer_rev	Marker	 Biosample	Replicate	Run	Fastq_fwd	Fastq_rev
    cgatcgtcatcacg	TCCACTAATCACAARGATATTGGTAC	cgcgatctgtagag	WACTAATCAATTWCCAAATCCTCC	MFZR	14Mon01	repl2	prerun	MFZR2_S5_L001_R1_001.fastq	MFZR2_S5_L001_R2_001.fastq

Merge FASTQ files
-------------------------------------------------------------------------

In addition to */path/to/my/fastqdir* and *fastqinfo.tsv*, we need

- Output WopMars DB file
- Output TSV file with FASTA file sample information
- Output directory to write the merged FASTA files

.. code-block:: bash

    vtam --wopdb wopdb.sqlite --fastqinfo fastqinfo.tsv --fastainfo fastainfo.tsv --fastqdir /path/to/my/fastqdir --fastadir /path/to/my/fastadir


Check that the *fastainfo.tsv* file is there.

.. code-block:: bash

    $ head -n1 fastainfo.tsv
    cgatcgtcatcacg	TCCACTAATCACAARGATATTGGTAC	cgcgatctgtagag	WACTAATCAATTWCCAAATCCTCC	MFZR	14Mon01	repl2	prerun	prerun_MFZR_repl2.fasta

and also the expected *FASTA* files

.. code-block:: bash

    $ ls  /path/to/my/fastadir
    prerun_MFZR_repl2.fasta  prerun_MFZR_repl3.fasta  prerun_ZFZR_repl1.fasta  prerun_ZFZR_repl2.fasta  prerun_ZFZR_repl3.fasta


CONTINUEEE

Run the *SampleInformation* rule
-------------------------------------------------------------------------

If the *Merge* rule was fine, then you can run until the *SampleInformation* step

.. code-block:: bash

    wopmars -w Wopfile.yml -D "sqlite:///db.sqlite" -v -p -t SampleInformation

Check the *db.sqlite* file using *sqlitebrowser*. All table except *Variant* should be filled.


Run the *VariantReadCount* rule
-------------------------------------------------------------------------

We can continue with the *VariantReadCount* rule. This rule is quite long, so be patient or test it first with smaller datasets.

.. code-block:: bash

    wopmars -w Wopfile.yml -D "sqlite:///db.sqlite" -v -p -t VariantReadCount

This rule fills in the *Variant* table with the number of reads per variant and marker across all samples.

You can have more count detail in each category in the tmp files listed in the *sortreads_samplecount.tsv* file


Run the *FilterMinReplicateNumber* rule
-------------------------------------------------------------------------

We can continue with the *FilterMinReplicateNumber* rule. This rule is quite long, so be patient or test it first with smaller datasets.

.. code-block:: bash

    wopmars -w Wopfile.yml -D "sqlite:///db.sqlite" -v -p -t FilterMinReplicateNumber

This rule fills in the *Variant* table with the number of reads per variant and marker across all samples.

You can have more count detail in each category in the tmp files listed in the *sortreads_samplecount.tsv* file


Run the *TaxAssign* rule
-------------------------------------------------------------------------

We can continue with the *TaxAssign* rule. This rule is long, so be patient or test it with smaller datasets.

.. code-block:: bash

    wopmars -w Wopfile.yml -D "sqlite:///db.sqlite" -v -p -t TaxAssign

This rule determine the LTG (Lower taxonomic group) for each variant

