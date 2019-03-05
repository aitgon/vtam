Developer - FilterMinReplicateNumber wrapper
=================================================

Input
---------------------

Output
---------------------

The FilterMinReplicateNumber wrapper outputs three files:

- filtered_dataframe_path.tsv
- MFZR_variant_info.tsv
- MFZR_variant.fasta

filtered_dataframe_path.tsv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The path of this file is defined in the Wopfile *FilterMinReplicateNumber* rule here:

.. code-block:: bash

    output:
        file:
            filtered_dataframe_path: data/filtered_dataframe_path.tsv

This file has three columns with 

- Marker_name
- Path to a TSV file with variants of a marker that passed the filters and additional information from the FilterMinReplicateNumber wrapper
- Path to a FASTA file with sequences of those variants

Example of *filtered_dataframe_path.tsv*

.. code-block:: bash

    MFZR	data/output/FilterMinReplicateNumber/MFZR_variant_info.tsv	data/output/FilterMinReplicateNumber/MFZR_variant.fasta
    ZFZR	data/output/FilterMinReplicateNumber/ZFZR_variant_info.tsv	data/output/FilterMinReplicateNumber/ZFZR_variant.fasta

The directory path of the TSV and FASTA files of the markers is defined in the Wopfile *FilterMinReplicateNumber* rule here:

.. code-block:: bash

    params:
        filter_output_dir: data/output/FilterMinReplicateNumber

First two lines of *MFZR_variant_info.tsv* for marker *MFZR*

.. code-block:: bash

	    variant_seq	replicate	biosample	sample_replicate	count	is_borderline	is_pseudogene_indel	is_pseudogene_codon_stop	read_average
    0	TTTATACTTCATTTTTGGGGCTTGATCTGGTATAGTAGGGACATCTCTTAGTCTACTAATTCGAGCTGAATTAGGACAACCAGGATCCCTTATTGGAGACGACCAAATTTACAATGTAATTGTCACAGCCCATGCCTTTATTATAATTTTCTTCATGGTTATGCCCATTATAATT	repl2	14Cro11	14Cro11_repl2	26	False	False	False	12.0

First two lines of *MFZR_variant.fasta* for marker *MFZR*

.. code-block:: bash

    >TTTATATTTTATTTTTGGTGTTTGAGCCGGAATAATTGGCTTAAGAATAAGCCTGCTGATCCGTTTAGAGCTTGGGGTTTTAGGGCCTTTTCTTGGAGACGAGCATTTGTATAACGTTATTGTTACTGCCCATGCTTTTGTTATAATTTTCTTTATAGTTATACCAATTTCTATA
    TTTATATTTTATTTTTGGTGTTTGAGCCGGAATAATTGGCTTAAGAATAAGCCTGCTGATCCGTTTAGAGCTTGGGGTTTTAGGGCCTTTTCTTGGAGACGAGCATTTGTATAACGTTATTGTTACTGCCCATGCTTTTGTTATAATTTTCTTTATAGTTATACCAATTTCTATA


