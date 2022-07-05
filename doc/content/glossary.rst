Glossary
=======================================


Terms are defined as they are used in VTAM. They might not hold or be sufficiently precise in another context outside VTAM. __


.. _ASV_glossary:

ASV or variant
------------------------------------------------

Amplicon Sequence Variant: Unique amplicon sequence (Callahan et al., 2017). Identical reads are pooled into a variant. Variants are characterized by the number of reads in each replicate of each sample, which we also call “sample-replicate”.


.. _ASVtable_glossary:

ASV table
------------------------------------------------

Representation of presence of each of the variants in each sample; Variants are in lines, samples are in columns, read numbers or presence absence are in cells.


.. _BLASThit_glossary:

BLAST hit
------------------------------------------------

A sequence from the BLAST database that has significant similarity to the query sequence (variant).

.. _borderline_glossary:

Chimera borderline
------------------------------------------------

When the chimera formation happens near the extremity of the parental sequences, the resulting chimera is very similar to one of the parental sequences. These chimeras are difficult to tell apart from real variants.


.. _coverage_glossary:

Coverage
------------------------------------------------

[in BLAST] the percentage of the length of the query sequence that is covered by the BLAST alignment.


.. _demultiplexing_glossary:

Demultiplexing
------------------------------------------------

Sorting reads to sample-replicates according to the presence of primers and tags at their extremities.


.. _dereplication_glossary:

Dereplication
------------------------------------------------

Identical reads are pooled into a variant, and the read count is kept as an information.


.. _flag_glossary:

Flag
------------------------------------------------

If variants or occurrences are flagged, they remain in the dataset after the corresponding filtering step, but they will be marked (flagged) in the ASV table.


.. _locus_glossary:

Locus/Gene
------------------------------------------------

genomic region (COI, LSU, MatK, RBCL...) 


.. _LTG_glossary:

Lowest Taxonomic Group (LTG)
------------------------------------------------

The taxonomic group of the highest resolution (species is high resolution, phylum is low), that contains all or a given % of the sequences.


.. _marker_glossary:

Marker
------------------------------------------------

A region amplified by one primer pair.

.. _merge_glossary:

Merge
------------------------------------------------

Assemble each forward and reverse reads (read pair) to a single sequence.


.. _mock_glossary:

Mock sample
------------------------------------------------

Sample with known DNA composition.


.. _modulo_glossary:

Modulo
------------------------------------------------

Remainder after a division of one number by another. (e.g. the modulo 3 of 4 is 1)


.. _occurrence_glossary:

Occurrence
------------------------------------------------

Presence of a variant in a sample or sample-replicate

.. _delete_glossary:

Unexpected or ‘delete’ occurrence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A variant in a sample that is known to be erroneous. It can be a variant in a negative control, an unexpected variant in a mock sample, or variant identified from a clearly different habitat than that of the sample.

.. _keep_glossary:

Expected or ‘keep’ occurrence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A variant that should be present in the given sample after filtering. These are expected variants in the mock samples.


.. _OLSP_glossary:

OLSP
------------------------------------------------

One-Locus-Several-Primers Strategy of using more than one primer pairs that amplify the same locus (slight variation in the position of the annealing sites) in order to increase the taxonomic coverage (Corse et al., 2019).


.. _renkonen_glossary:

Renkonen distance
------------------------------------------------

1 - sum(p1i, p2i), where p1i is the frequency of variant i in sample-replicate1 (Renkonen, 1938)


.. _replicate_glossary:

Replicate 
------------------------------------------------

or Replicate series : Pool of a single replicate from each sample of a run.


.. _run_glossary:

Run
------------------------------------------------

A pool of samples and the associated positive (mock) and negative controls. Ideally, they are obtained in the same sequencing run. 


.. _sample_glossary:

Sample
------------------------------------------------

DNA extraction from a given environment/individual.

.. _sample-replicate_glossary:

Sample-Replicate
------------------------------------------------

Technical replicate of the same sample. e.g different PCR on the same DNA extraction.


.. _tag-jump_glossary:

Tag-jump
------------------------------------------------

Generation of artefactual sequences in which amplicons carry different tags than originally applied (Schnell et al., 2015)


.. _tag_glossary:

Tag
------------------------------------------------

Short DNA sequences present at one or both extremities of the amplified DNA fragment. A tag or the combination of forward and reverse tags determine the sample-replicate where the read comes from.


.. _trimming_glossary:

Trimming
------------------------------------------------

Removing part of the extremities of a sequence (e.g. trim the tags/adapters/primers from a read to obtain the biological sequence)


.. _TSV_glossary:

TSV
------------------------------------------------

A text file format with tab-separated values.


