Input/Output files
============================

Help on the usage and complete list of I/O arguments of each command can be obtained using the command line help

.. code-block:: bash

    vtam COMMAND --help
    i.e.
    vtam filter --help

Here we detail the content of the I/O files

.. _params_io:

params
--------------------------------

Input of most commands. YML file with :ref:`numerical parameters <numerical_parameter_file_reference>`. Can be omitted if all parameters are by default.
Simple text file with a “parameter name: parameter value” format. One parameter per line
e.g. 

.. code-block:: bash

    lfn_variant_cutoff: 0.001
    lfn_sample_replicate_cutoff: 0.003
    lfn_read_count_cutoff: 70
    pcr_error_var_prop: 0.05

.. _fastqinfo_io:

fastqinfo
--------------------------------

Input of :ref:`merge <merge_reference>`. TSV file with the following columns:

    - TagFwd: Sequence of the tag on the forward primer (5’=>3’)
    - PrimerFwd: Sequence of the forward primer (5’=>3’)       
    - TagRev: Sequence of the tag on the reverse primer (5’=>3’)
    - PrimerRev: Sequence of the reverse primer (5’=>3’) 
    - Marker: Name of the marker (e.g. MFZR)
    - Sample: Name of the sample
    - Replicate: ID of the replicate 
    - Run: Name of the sequencing run
    - FastqFwd: Name of the forward fastq file
    - FastqRev: Name of the reverse fastq file

.. _fastainfo_io:

fastainfo
--------------------------------

Output of :ref:`merge <merge_reference>`, input of :ref:`sortreads <sortreads_reference>`. TSV file with the following columns:

    - run: Name of the sequencing run
    - marker: Name of the marker (e.g. MFZR) 
    - sample: Name of the sample
    - replicate: ID of the replicate
    - tagfwd: Sequence of the tag on the forward primer (5’=>3’)
    - primerfwd: Sequence of the forward primer (5’=>3’)
    - tagrev: Sequence of the tag on the revrese primer (5’=>3’) 
    - primerrev: Sequence of the reverse primer (5’=>3’) 
    - mergedfasta: name of the fasta file with merged sequences

.. _sortedinfo_io:

sortedinfo
--------------------------------

Output of :ref:`sortreads <sortreads_reference>`, input of :ref:`filter <filter_reference>` and :ref:`optimize <optimize_reference>`. TSV file with the following columns:

    - run: Name of the sequencing run
    - marker: Name of the marker (e.g. MFZR) 
    - sample: Name of the sample
    - replicate: ID of the replicate
    - sortedfasta: name of the fasta file containing merged, demultiplexed, trimmed sequeces

.. _db_io:

db
--------------------------------

I/O of :ref:`filter <filter_reference>`, :ref:`taxassign <taxassign_reference>`. Input of :ref:`optimize <optimize_reference>`, :ref:`pool <pool_reference>`. 
Sqlite database containing variants, samples, replicates, read counts, information on filtering steps, taxonomic assignations.

.. _asvtable_io:

asvtable
--------------------------------

Output of :ref:`filter <filter_reference>` or :ref:`pool <pool_reference>`, input of :ref:`taxassign <taxassign_reference>`. 
TSV file with the variants (in lines) that passed all filtering steps, samples (in columns), presence-absence (output of pool) or read counts (output of filter) in cells and additional columns:

    - run: Name of the sequencing run
    - marker: Name of the marker (e.g. MFZR)   
    - variant: Variant ID
    - pooled_variants (only in output of pool): IDs of variants pooled since identical in their overlapping regions
    - sequence_length length of the variant
    - read_count: Total number of reads of the variants in the samples listed in the table
    - [one column per sample] :  presence-absence (output of pool) or read counts (output of filter)
    - clusterid: ID of the centroïd of the cluster (0.97 clustering of all variants of the asv table)
    - clustersize: Number of variants in the cluster
    - chimera_borderline (only in output of filter): Potential chimeras (very similar to one of the parental sequence)
    - [keep_mockXX; One column per mock sample, if known_occurrences option is used]: 1 if variant is expected in the mock sample, 0 otherwise
    - pooled_sequences (only in output of pool): Sequences of pooled_variants
    - sequence: Sequence of the variant

.. _known_occurrences_io:

:ref:`known_occurrences <optimize_reference>`
----------------------------------------------------------------

Input of :ref:`filter <filter_reference>` and :ref:`optimize <optimize_reference>`. TSV file with expected occurrences (keep) and known false positives (delete). 

    - Marker: Name of the marker (e.g. MFZR) 
    - Run: Name of the sequencing run
    - Sample: Name of the sample
    - Mock: 1 if sample is a mock, 0 otherwise
    - Variant: Varinat ID (can be empty)
    - Action: keep (occurrences that should be kept after filtering) or delete (clear false positives)
    - Sequence: Sequence of the variant
    - Tax_name: optional, not used by optimize

.. _optimize_lfn_sample_replicate_io:

:ref:`optimize_lfn_sample_replicate.tsv <OptimizeLFNsampleReplicate_reference>`
------------------------------------------------------------------------------------------------

Output of :ref:`optimize <optimize_reference>`. TSV file with the following columns:

    - run: Name of the sequencing run
    - marker: Name of the marker (e.g. MFZR) 
    - sample: Name of the sample
    - replicate: ID of the replicate
    - variant: Variant ID
    - N_ijk: Number of reads of variant i, in sample j and replicate k
    - N_jk: Number of reads in sample j and replicate k (all variants)
    - N_ijk/N_jk
    - round_down: Rounded value of N_ijk/N_jk
    - sequence: Variant sequence

.. _optimize_lfn_read_count_and_lfn_variant_io:

:ref:`optimize_lfn_read_count_and_lfn_variant.tsv OR optimize_lfn_read_count_and_lfn_variant_replicate.tsv <OptimizeLFNReadCountAndLFNvariant_reference>`
----------------------------------------------------------------------------------------------------------------------------------------------------------------

Output of :ref:`optimize <optimize_reference>`. TSV file with the following columns:

    - occurrence_nb_keep: Number of keep occurrence left after filtering with lfn_nijk_cutoff and lfn_variant_cutoff values
    - occurrence_nb_delete: Number of delete occurrence left after filtering with lfn_nijk_cutoff and lfn_variant_cutoff values
    - lfn_nijk_cutoff: lfn_read_count_cutoff
    - lfn_variant_cutoff or lfn_variant_replicate_cutoff
    - run: Name of the sequencing run
    - marker: Name of the marker (e.g. MFZR) 

.. _optimize_lfn_variant_specific_io:

:ref:`optimize_lfn_variant_specific.tsv OR optimize_lfn_variant_replicate_specific.tsv <OptimizeLFNReadCountAndLFNvariant_reference>`
-------------------------------------------------------------------------------------------------------------------------------------

Output of :ref:`optimize <optimize_reference>`. TSV file with the following columns:

    - run: Name of the sequencing run
    - marker: Name of the marker (e.g. MFZR)  
    - variant: Variant ID
    - replicate: (if optimize_lfn_variant_replicate_specific.tsv) ID of the replicate
    - action: Type d’occurrece (delete/keep) 
    - read_count_max: Max of N_ijk for a given i
    - N_i (optimize_lfn_variant_specific.tsv) : Number of reads of variant i 
    - N_ik (optimize_lfn_variant_replicate_specific.tsv): Number of reads of variant i in replicate k
    - lfn_variant_cutoff: read_count_max/N_i or read_count_max/N_ik 
    - sequence: Variant sequence

.. _optimize_pcr_error_io:

:ref:`optimize_pcr_error.tsv <OptimizePCRError_reference>`
------------------------------------------------------------------------------------------------

Output of :ref:`optimize <optimize_reference>`. TSV file with the following columns:

    - run: Name of the sequencing run
    - marker: Name of the marker (e.g. MFZR) 
    - sample: Name of the sample
    - variant_expected: ID of a keep variant
    - N_ij_expected: Number of reads of the expected variant in the sample (all replicates)
    - variant_unexpected: ID of an unexpected variants with one mismatch to the keep variant
    - N_ij_unexpected: Number of reads of the unexpected variant in the sample (all replicates)
    - N_ij_unexpected_to_expected_ratio: N_ij_unexpected/N_ij_expected
    - sequence_expected: Sequence of the expected variant
    - sequence_unexpected: Sequence of the unexpected variant

.. _output_io:

output (taxassign)
--------------------------------

Output of :ref:`taxassign <taxassign_reference>`
The input asvtable completed with the following columns:

    - ltg_tax_id: TaxID of the LTG (Lowest Taxonomic Group)
    - ltg_tax_name    ltg_rank: Name of the LTG
    - identity: Percentage of identity used to determine the LTG
    - blast_db: Name of the taxonomic BLAST database files (without extensions)
    - phylum: Phylum of LTG
    - class: class  of LTG
    - order: order  of LTG
    - family: family  of LTG
    - genus: genus of LTG
    - species: species  of LTG

.. _taxonomy_io:

taxonomy
--------------------------------

Output of :ref:`taxonomy <taxonomy_reference>`, input of :ref:`taxassign <taxassign_reference>`. TSV file with information of all taxa in the reference (BLAST) database.

    - tax_id: Taxonomic identifier of the taxon
    - parent_tax_id: Taxonomic identifier of the direct parent of the taxon
    - rank: Taxonomic rank of the taxon (e.g. class, species, no rank)
    - name_txt: Name of the taxon
    - old_tax_id: TaxID of taxa merged to taxon (not valid any more)

.. _runmarker_io:

runmarker
--------------------------------

Input of :ref:`pool <pool_reference>`. TSV file with the list of all run-marker combinations to be pooled.

    - run: Name of the sequencing run
    - marker: Name of the marker (e.g. MFZR) 
