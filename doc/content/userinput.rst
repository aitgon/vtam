User Input
===================================

The user needs to provide two CSV files

- FASTQ files
- A CSV file with information about the FASTQ files
- A CSV file with information about the samples

CSV file with information about the FASTQ files
----------------------------------------------------

It takes this format:

.. code-block:: bash

    fastq_fw;fastq_rev;output_fas
    test1.fastq;test2.fastq;prerun_test.fasta

CSV file with information about the samples
----------------------------------------------------

It takes this format

.. code-block:: bash

    Tag Forward,Primer Forward,Tag Reverse,Primer Reverse,Marker name ,Sample,Replicate,Filename,Run
    cgatcgtcatcacg,TCCACTAATCACAARGATATTGGTAC,cacgatttgtagag,WACTAATCAATTWCCAAATCCTCC,MFZR,14Mon01,repl1,data/fastq_merged_fasta/prerun_MFZR_repl1.fasta,prerun
    tcgatcacgatgt,TCCACTAATCACAARGATATTGGTAC,cgcgtctgtagag,WACTAATCAATTWCCAAATCCTCC,MFZR,14Mon02,repl1,data/fastq_merged_fasta/prerun_MFZR_repl1.fasta,prerun

The "avant dernier" column in this file must correspond to the some output of the merge step


