Tutorial
============

To setup a new Wopmetabarcodign workflow, follow these steps:

1. Download Wopfile
2. Check FASTQ location
3. Specify FASTQ pair input and FASTA output information (*sample2fastq.csv*)
4. Specify sample information with FASTA paths (*sample2fasta.csv*)
5. Update *sample2fasta.csv* paths in Wopfile
6. Update workflow parameters in Wopfile (Optional)
7. Run wopmars

Download Wopfile or Workflow definition file
---------------------------------------------------

This is a full :download:`Wopfile example </data/Wopfile.yml>`. Here there are the first lines:

.. code-block:: bash

    rule Merge:
      tool: wopmetabarcoding.wrapper.Merge
      input:
          file:
              sample2fastq: data/sample2fastq_test.csv
      output:
          file:
              sample2fasta: data/sample2fasta.csv
      params:
          fastq_directory: data/input/fastq
    ...

Check FASTQ file location
---------------------------------------------------------------------------------

In my case, FASTQ files are here:

.. code-block:: bash

    $ ls /home/gonzalez/Data/2017_meglecz_metabarcoding/data/nuxeo/input/input_fastq
    MFZR1_S4_L001_R1_001.fastq  MFZR2_S5_L001_R2_001.fastq  ZFZR1_S1_L001_R1_001.fastq  ZFZR2_S2_L001_R2_001.fastq
    MFZR1_S4_L001_R2_001.fastq  MFZR3_S6_L001_R1_001.fastq  ZFZR1_S1_L001_R2_001.fastq  ZFZR3_S3_L001_R1_001.fastq
    MFZR2_S5_L001_R1_001.fastq  MFZR3_S6_L001_R2_001.fastq  ZFZR2_S2_L001_R1_001.fastq  ZFZR3_S3_L001_R2_001.fastq


Set FASTQ input and FASTA output directories in Wopfile
---------------------------------------------------------------------

Based on my FASTQ file location and desired FASTA output location, I set in the *Merge* rule of the *Wopfile*:

.. code-block:: bash

    params:
        fastq_directory: /home/gonzalez/Data/2017_meglecz_metabarcoding/data/nuxeo/input/input_fastq
        fasta_output_dir: /home/gonzalez/Data/2017_meglecz_metabarcoding/repositories/wopmetabarcoding-appli/out/fasta


Define sample and FASTQ file pair information and set path in Wopfile
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

The first two lines of my *sample2fastq.tsv* look like this:

.. code-block:: bash

    Tag_fwd	Primer_fwd	Tag_rev	Primer_rev	Marker	 Biosample	Replicate	Run	Fastq_fwd	Fastq_rev
    cgatcgtcatcacg	TCCACTAATCACAARGATATTGGTAC	cgcgatctgtagag	WACTAATCAATTWCCAAATCCTCC	MFZR	14Mon01	repl2	prerun	MFZR2_S5_L001_R1_001.fastq	MFZR2_S5_L001_R2_001.fastq

Set path to *sample2fastq.tsv* file in Wopfile. I set in the *Merge* rule of the *Wopfile*:

.. code-block:: bash

  input:
      file:
          sample2fastq: /home/gonzalez/Data/2017_meglecz_metabarcoding/repositories/wopmetabarcoding-appli/sample2fastq.tsv

Run the *Merge* rule
-------------------------------------------------------------------------

At this point, it is already possible to test the *Merge* rule. First try a dry-run using the Wopmars '-n' option

.. code-block:: bash

    wopmars -w Wopfile.yml -D "sqlite:///db.sqlite" -v -p -t Merge -n


And then a true run without the '-n' option

.. code-block:: bash

    wopmars -w Wopfile.yml -D "sqlite:///db.sqlite" -v -p -t Merge

Check that the *sample2fasta.tsv* file is there.

.. code-block:: bash

    $ head -n1 out/sample2fasta.tsv
    cgatcgtcatcacg	TCCACTAATCACAARGATATTGGTAC	cgcgatctgtagag	WACTAATCAATTWCCAAATCCTCC	MFZR	14Mon01	repl2	prerun	prerun_MFZR_repl2.fasta

and also the expected *FASTA* files

.. code-block:: bash

    $ ls  out/fasta/
    prerun_MFZR_repl2.fasta  prerun_MFZR_repl3.fasta  prerun_ZFZR_repl1.fasta  prerun_ZFZR_repl2.fasta  prerun_ZFZR_repl3.fasta


Run the *SampleInformation* rule
-------------------------------------------------------------------------

If the *Merge* rule was fine, then you can run until the *SampleInformation* step

.. code-block:: bash

    wopmars -w Wopfile.yml -D "sqlite:///db.sqlite" -v -p -t SampleInformation

Check the *db.sqlite* file using *sqlitebrowser*. All table except *Variant* should be filled.


Run the *SortReads* rule
-------------------------------------------------------------------------

We can continue with the *SortReads* rule. This rule is quite long, so be patient or test it first with smaller datasets.

.. code-block:: bash

    wopmars -w Wopfile.yml -D "sqlite:///db.sqlite" -v -p -t SortReads

This rule fills in the *Variant* table with the number of reads per variant and marker across all samples.

You can have more count detail in each category in the tmp files listed in the *sortreads_samplecount.tsv* file


Run the *Filter* rule
-------------------------------------------------------------------------

We can continue with the *Filter* rule. This rule is quite long, so be patient or test it first with smaller datasets.

.. code-block:: bash

    wopmars -w Wopfile.yml -D "sqlite:///db.sqlite" -v -p -t Filter

This rule fills in the *Variant* table with the number of reads per variant and marker across all samples.

You can have more count detail in each category in the tmp files listed in the *sortreads_samplecount.tsv* file


Run the *Taxassign* rule
-------------------------------------------------------------------------

We can continue with the *Taxassign* rule. This rule is long, so be patient or test it with smaller datasets.

.. code-block:: bash

    wopmars -w Wopfile.yml -D "sqlite:///db.sqlite" -v -p -t Taxassign

This rule determine the LTG (Lower taxonomic group) for each variant

