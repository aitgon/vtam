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

The "avant dernier" column in this file must correspond to the some output of the merge step.

CSV file with cutoffs values for filters
----------------------------------------------------

.. code-block:: bash

    sequence	value
    CCTTTATCTAGTATTCGGTGCTTGGGCTGGGATAGTTGGAACAGCCCTTAGCTTACTAATCCGTGCAGAGCTTAGCCAACCTGGCGCCCTGCTCGGTGACGACCAAGTTTACAACGTGATCGTAACAGCTCATGCTTTCGTAATAATCTTCTTTATAGTAATGCCAATTATGATT	0.002
    TTTATATTTTATTTTTGGAGCCTGAGCTGGAATAGTAGGTACTTCCCTTAGTATACTTATTCGAGCCGAATTAGGACACCCAGGATCTCTAATTGGAGACGACCAAATTTATAATGTAATTGTTACTGCTCATGCTTTTGTAATAATTTTTTTTATAGTTATACCAATTATAATT	0.003

CSV file with genetic codes
---------------------------------------------------

.. code-block:: bash

    TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG;Base1
    TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG;Base2
    TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG;Base3

    FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG;1;The Standard Code (transl_table=1)
    FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG;2;The Vertebrate Mitochondrial Code (transl_table=2)

CSV file with taxonomic association level
---------------------------------------------------

.. code-block:: bash

    100.0	species	subspecies	1
    97.0	genus	species	1
    95.0	family	species	3
    90.0	order	family	3
    85.0	order	order	3
    80.0	class	order	5

Reference fasta database for taxassign vsearch
---------------------------------------------------

.. code-block:: bash

    >5244419 name=Echinorhynchida tax_id=57283 rank=order parent_taxid=45080
    TTAATGTATGTGTTAGTGGGAGTTTGAGGGGGGTTGATAGGATTCTCAATAAGGATGTTAATTCGTTTAGAGCTAGGAAGTGGGGGCATTTGAATGGGTAGG

