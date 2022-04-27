Reference
===================================

.. _numerical_parameter_file_reference:

The numerical parameter file
--------------------------------

For each vtam command the **params** argument takes an optional YML file that sets parameter values. Omitting the **params** argument will prompt vtam to use default values. The full list of parameters and the default values is shown here:

.. code-block:: bash

	################################################################################
	# Parameters of the "merge" command
	# These parameters are used by the vsearch --fastq_mergepairs tool 
	# that underlies the "vtam merge" command
	# For a description of these parameters run_name "vsearch --help"
	fastq_ascii: 33
	fastq_maxee: 1
	fastq_maxmergelen: 500
	fastq_maxns: 0
	fastq_minlen: 50
	fastq_minmergelen: 100
	fastq_minovlen: 50
	fastq_truncqual: 10
	 
	################################################################################
	# Parameters of the "sortreads" command
	# These parameters correspond to the corresponding parametes by cutadapt 
	# that underlies the "vtam sortreads" command
	# For a description of these parameters run_name "cutadapt --help"
	cutadapt_error_rate: 0.1 # -e in cutadapt
	cutadapt_minimum_length: 50 # -m in cutadapt
	cutadapt_maximum_length: 500 # -M in cutadapt
	 
	################################################################################
	# Parameters of the "filter" command
	# This parameter sets the minimal number of reads of a variant in the whole run
	global_read_count_cutoff: 2
	 
	################################################################################
	# Parameters of the "FilterLFN" filter in the "filter" command
	# These parameters set the cutoffs for the low frequency noise (LFN) filters
	# N_ijk is the number of reads of varinat i, sample j and replicate k
	# Occurrence is deleted if N_ijk/N_i < lfn_variant_cutoff
	lfn_variant_cutoff: 0.001
	# Occurrence is deleted if N_ijk/N_ik < lfn_variant_replicate_cutoff
	# This parameter is used if the --lfn_variant_replicate option is set in 
	# "vtam filter" or "vtam optimize"
	lfn_variant_replicate_cutoff: 0.001
	# Occurrence is deleted if N_ijk/N_jk < lfn_ sample lfn_sample_replicate_cutoff
	lfn_sample_replicate_cutoff: 0.001
	# Occurrence is deleted if N_ijk < lfn_ lfn_read_count_cutoff
	lfn_read_count_cutoff: 10
	 
	################################################################################
	# Parameters of the "FilterMinReplicateNumber" filter in the "filter" command
	# Occurrences of a variant in a given sample are retained only if it is 
	# present in at least min_replicate_number replicates of the sample
	min_replicate_number: 2
	 
	################################################################################
	# Parameter of the "FilterPCRerror" filter in the "filter" command
	# A given variant 1 is eliminated if N_1j/N_2j < pcr_error_var_prop, where 
	# variant 2 is identical to variant 1 except a single mismatch
	pcr_error_var_prop: 0.1
	 
	################################################################################
	# Parameter of the "FilterChimera" filter in the "filter" command
	# This parameter corresponds to the abskew parameter in the vsearch 
	# --uchime3_denovo tool that underlies the vtam FilterChimera
	# For a description of this parameter run_name "vsearch --help"
	uchime3_denovo_abskew: 16.0
	 
	################################################################################
	# Parameter of the "FilterRenkonen" filter in the "filter" command
	# Quantile renkonen distance to drop more extreme values
	# For. a 0.9 value will set the 9th decile of all renkonen distances as cutoff
	renkonen_distance_quantile: 0.9
	 
	################################################################################
	# Parameter of the "FilterIndel" filter in the "filter" command
	# If 1, skips this filter for non-coding markers
	skip_filter_indel: 0
	 
	################################################################################
	# Parameter of the "FilterCondonStop" filter in the "filter" command
	# If 1, skips this filter for non-coding markers
	skip_filter_codon_stop: 0
	# Translation table number from NCBI [ link]
	# Default NCBI translation table 5: stops: ['TAA', 'UAA', 'TAG', 'UAG']
	genetic_code: 5
	 
	################################################################################
	# Parameter of the "MakeAsvTable" filter in the "filter" command
	# Cluster identity value to clusterize sequences
	cluster_identity: 0.97
	 
	################################################################################
	# Parameters of the "taxassign" command
	# Blast parameter for the minimum query coverage
	qcov_hsp_perc: 80
	# The LTG must include include_prop percent of the hits
	include_prop: 90
	# Minimal number of taxa among the hits to assign LTG when %identity 
	# is below ltg_rule_threshold
	min_number_of_taxa: 3
	ltg_rule_threshold: 97

.. _merge_reference:

The command merge
--------------------------------

VTAM can start from FASTQ files of paired end metabarcoding data. This command merges paired end sequences using the **fastq_mergepairs** command of `vsearch <https://github.com/torognes/vsearch>`_. 

A quick introduction to the **vtam merge** command is given in the tutorial :ref:`merge_tutorial`. 

The command line arguments of the **merge** command can be obtained with 

.. code-block:: bash

    vtam merge --help

The most important arguments are these inputs and outputs:

Inputs:
    - :ref:`fastqinfo <fastqinfo_io>`: TSV file with files to be merged. These files can be compressed in .gz or .bz2.
    - **fastqdir**: Path to the directory containing the fastq files


Outputs:
    - :ref:`fastainfo <fastainfo_io>`: TSV file created by vtam merge. It contains all the info of :ref:`fastqinfo <fastqinfo_io>` completed by the names of the merged fasta files.
    - **fastadir**: Directory to keep the output merged fasta files.

The list of numerical parameters can be found in the `The numerical parameter file`_ section.

.. _sortreads_reference:

The command sortreads
--------------------------------

Typically, the sequencing reads contain primers and tags, and this command uses them to attribute each read to a run-marker-sample-replicate. 

A quick introduction to the **sortreads** command is given in the tutorial :ref:`sortreads_tutorial`.

The arguments of the **sortreads** command can be obtained with 

.. code-block:: bash

    vtam sortreads --help

The most important arguments are these inputs and outputs:

Inputs:
    - :ref:`fastainfo <fastainfo_io>`: TSV file created by **merge** (allows gzip and bzip2 compressed files). It contains all the info of :ref:`fastqinfo <fastqinfo_io>` completed by the names of the merged fasta files.
    - **fastadir**: Directory containing the merged fasta files.

Outputs:
    - **sorteddir**: Directory to keep the output :ref:`demultiplexed <demultiplexing_glossary>` fasta files. In this folder, there is also a :ref:`sortedinfo <sortedinfo_io>` file with the information about each sorted fasta file.

The list of numerical parameters can be found in the `The numerical parameter file`_ section. 

These are the different actions of the **sortreads** command:
    - Sort reads to sample-replicates according to the presence of tags using cutadapt. Exact matches are imposed between reads and tags, and the minimum overlap is the length of the tag.
    - Trim reads from primers using cutadapt. Mismatches are allowed (**cutadapt_error_rate**), but no indels. The minimum overlap between the primer and the read is the length of the primer.
    - The trimmed sequences are kept if their length is between **cutadapt_minimum_length** and **cutadapt_maximum_length**.

The optionnal argument **no_reverse** can be used if the dataset does not contain any reverse sequence. It makes the run faster.

The optionnal argument **tag_to_end** can be used if the sequences' tags are located at the edges of the sequences.It make the run faster.

The optionnal argument **primer_to_end** can be used if the sequencess primers are located at the edges of the sequences (between the tags and the primers). It makes the run faster.



.. random_seq:

The command random_seq (OPTIONNAL)
--------------------------------

When working with large datasets VTAM can subselect a random set of sequences in order to run with less data to process to reduce the running time and the workoad on a users machine.

A quick introduction to the **vtam random_seq** command is given in the tutorial :ref:`random_seq_tutorial`. 

The command line arguments of the **random_seq** command can be obtained with 

.. code-block:: bash

    vtam random_seq --help

The most important arguments are these inputs and outputs:

Inputs:
    - :ref:`fastainfo <fastainfo_io>`: TSV file with files to take the sequences from
    - **fastadir**: Path to the directory containing the fasta files
    - **samplesize**: number of sequences to be randomly selected


Outputs:
    - :ref:`random_seqinfo <random_seqinfo_io>`: TSV file created by vtam random_seq. It contains all the info of :ref:`fastainfo <fastqinfo_io>` completed by the names of the fasta files containing the randomly selected seqquences.
    - **random_seqdir**: directory with randomly selected sequences in FASTA format



.. _filter_reference:

The command filter
--------------------------------

This command chains several steps of the analyses. It starts from :ref:`demultiplexed <demultiplexing_glossary>`, 
:ref:`trimmed <trimming_glossary>` reads, fills the sqlite database with variants and read counts, runs several 
filtering steps and produces :ref:`ASV tables <ASVtable_glossary>`. Each run-marker combination is treated independently during all filtering steps, even if several run-marker combinations are included in the dataset. 

A quick introduction to the **filter** command is given in the tutorial :ref:`filter_tutorial`.

The arguments of the **filter** command can be obtained with 

.. code-block:: bash

    vtam filter --help


These arguments are the most important inputs and outputs.
    
Inputs:
    - :ref:`sortedinfo <sortedinfo_io>`: TSV with the information about each sorted read file (output of sortreads). The content of the file will define the dataset to be used for filtering.
    - **sorteddir**: Directory containing the demultiplexed fasta files.

Database:
    - :ref:`db <db_io>`: Name of an sqlite db. It is used to store and access information on runs, markers, samples, variants, filtering steps and taxonomic assignations. It is created by **filter** if it does not exist and it is completed if it exists already. 
Outputs:
    - :ref:`asvtable <asvtable_io>`: TSV file containing variants that passed all the filters, together with read count in the different samples.

The list of numerical parameters can be found in the `The numerical parameter file`_ section.

In a typical run, **filter** should be run twice:
Once with light filtering parameters (default) to do a prefiltering step and produce a user-friendly :ref:`ASV table <ASVtable_glossary>`. 
Based on this ASV table users will identify the clearly :ref:`expected <keep_glossary>` and :ref:`unexpected occurrences <delete_glossary>` (e.g. in negative controls and mock samples). 
Using these occurrences the **optimize** command will suggest a parameter setting that keeps all expected occurrences and eliminates most unexpected ones. Then **filter** should be run again with the optimized parameter settings.

.. _make_known_occurrences_reference:

The command make_known_occurrences (OPTIONNAL)
--------------------------------

VTAM can create tsv file containing known and missing occurences base on an ASV table, a tsv file containing the samples types and a tsv file containing a mock composition
A quick introduction to the **vtam merge** command is given in the tutorial :ref:`merge_tutorial`. 

The command line arguments of the **make_known_occurrences** command can be obtained with 

.. code-block:: bash

    vtam make_known_occurrences --help

The most important arguments are these inputs and outputs:

Inputs:
    - :ref:`asvtable <ASVtable_io>`: TSV file containing variants in the last column (using ‘sequence’ as a header). Typically it is an ASV table. 
    - **sample_types**: Path to the tsv file containing the sample types
    - **mock_composition**: Path to the tsv file containing the mock composition


Outputs:
    - **known_occurrences**: path to the file containing the known occurrences.
    - **missing_occurrences**: path to the file containing the missing occurrences

.. _FilterLFN_reference:

FilterLFN
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This step intends to eliminate occurrences with low read counts that can be present due to light contaminations, sequencing errors and tag jumps.

Let *N_ijk* be the number of the reads of variant *i*, in sample *j* and replicate *k*.

Each :ref:`occurrence <occurrence_glossary>` can be characterized by the number of the reads of a variant in a given sample-replicate (*N_ijk*). Low read counts can be due to contamination or artefacts and therefore considered as Low Frequency Noise (LFN). 

The following LFN filters can be run on each occurrence, and the occurrence is retained only if it passes all of the activated filters:

    - LFN_sample_replicate filter: occurrence is deleted if *N_ijk*/*N_jk* < lfn_sample_replicate_cutoff
    - LFN_read_count filter: occurrence is deleted if *N_ijk* < lfn_read_count_cutoff
    - LFN_variant filter: occurrence is deleted if *N_ijk*/*N_i* < lfn_variant_cutoff
    - LFN_variant_replicate filter: occurrence is deleted if *N_ijk*/*N_ik* < lfn_variant_replicate_cutoff

**LFN_sample_replicate** and **LFN_read_count** filters intend to eliminate mainly sequencing or PCR artefacts.
**LFN_variant** and **LFN_variant_replicate filters** intend to eliminate occurrences that are present due to tag jump or slight inter sample contamination.

The **LFN_variant** and **LFN_variant_replicate** are two alternatives for the same idea and they are mutually exclusive. The **LFN_variant** mode is activated by default. Users can change to **LFN_variant_replicate** mode by using the **lfn_variant_replicate** flag in the command line. 
These filters eliminate occurrences that have low frequencies compared to the total number of reads of the variant (**LFN_variant**) or variant-replicate (**LFN_variant_replicate**) in the whole run-marker set. 

In a typical case, the **lfn_variant_cutoff** and **lfn_variant_replicate_cutoff** are the same for all variants for a given run-marker combination. This use should be preferred. However, occasionally, it can be justified to set individual (variant specific) thresholds to some of the variants. This is done using the **cutoff_specific** parameter that takes as a value a TSV file containing the variant specific threshold values. For variants not specified in the file, the value set by **lfn_variant_cutoff** or **lfn_variant_replicate_cutoff** is used.

.. _FilterMinReplicateNumber_reference:

FilterMinReplicateNumber
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**FilterMinReplicateNumber** is used to retain occurrences only if they are repeatable. 

Within each sample, occurrences of a variant are retained only if it is present in at least **min_replicate_number** replicates.

.. _FilterPCRerror_reference:

FilterPCRerror
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This step intends to eliminate occurrences due to PCR errors. These errors can be more frequent than sequencing errors.

Within each sample, this filter eliminates variants that have only one mismatch compared to another more frequent variant. The read count proportion of the two variants must be below **pcr_error_var_prop** in order to eliminate the least frequent variant.

.. _FilterChimera_reference:

FilterChimera
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Chimera filtering is done by the **uchime3_denovo** command integrated in `vsearch <https://github.com/torognes/vsearch>`_. Chimera checking is done sample by sample. Variants classed as chimeras are eliminated. Those classed as :ref:`borderline <borderline_glossary>` in at least one sample are flagged and this information will appear in the final ASV table.

.. _FilterRenkonen_reference:

The filter FilterRenkonen
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This step removes replicates that are very different to other replicates of the same sample and only makes sense if at least two replicates are used. The Renkonen distance takes into account the whole composition of the replicates.

Renkonen distances (:ref:`Renkonen, 1938 <Renkonen_1938_reflist>`) are calculated between each pair of replicates within each sample. Then a cutoff distance is set at a quantile of all renkonen distances of the dataset. This quantile can be given by the user through the **renkonen_distance_quantile** parameter and takes 0.9 (Ninth decile) as default value. Replicates above the quantile distance are removed.

.. _FilterIndel_reference:

FilterIndel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This filter makes only sense for coding sequences and can be skipped with the parameter **skip_filter_indel**.

It is based on the idea that sequences with indels of 3 nucleotides or its multiples are viable, but all others have frame-shift mutations and are unlikely to come from correct, functional sequences. Therefore, the :ref:`modulo <modulo_glossary>` 3 of the length of each variant is determined. The majority of the variants length will have the same modulo 3. All other variants are discarded.

.. _FilterCodonStop_reference:

FilterCodonStop
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This filter makes only sense for coding sequences and can be skipped with the parameter **skip_filter_codon_stop**. 

Given the appropriate `genetic code <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>`_, the presence of codon stops is checked in all reading frames of the direct strand. Variants that have a codon stop in all reading frames are discarded. 

When the Genetic code cannot be determined in advance (community includes species from groups using different genetic codes) we suggest using the Invertebrate Mithochondrial Code **genetic_table_number 5**, since its codons STOPs (TAA, TAG) are also codon STOP in almost all genetic codes.

.. _MakeAsvTable_reference:

MakeAsvTable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ASVs that passed all filters (Table *FilterCodonStop* in the SQLITE database) are written to an ASV table with variants in rows, samples in columns and read count per variants per samples are in cells. The read count per sample is the sum of read counts of the replicates. Replicates are not shown.

Additional columns include marker and run names, sequence length, total read counts for each variant, flags on :ref:`chimera borderline <borderline_glossary>` information on clustering and expected variants in mock samples (optional). 

**Clustering information**: 

When creating an ASV table, all variants are clustered using a user defined identity cutoff (**cluster_identity**; 0.97 by default). The **clusterid** is the centroid, and the number of variants are given in the **clustersize** column. This clustering is a simple help for the users to identify similar variants, but it is NOT part of the filtering process. Clusters of different ASV tables are not directly comparable. 

**Expected variants in mock samples**:

If mock samples have many expected variants or there are several different mock samples, it can be a bit difficult to identify all ‘keep’ and ‘delete’ :ref:`occurrences <occurrence_glossary>`. The **known_occurrences** option in filter command can help you in this task.

First, identify ‘keep’ occurrences in the mock samples, and make a *known_keep_occurences.tsv* file. Then run the filter command with the **known_occurrences** option. This will add an extra column foreach mock sample in the ASV table. These columns will contain 1 if the variant is expected in the mock sample. From this updated ASV table it will be easier to select ‘delete’ occurrences with the help of a tableur (LibreOffice, Excel). 

.. code-block:: bash

    vtam filter --db db.sqlite --sortedinfo sortedinfo.tsv --sorteddir sorted --asvtable asvtable_default.tsv -v --log vtam.log --known_occurrences known_keep_occurrences.tsv

The example of of the completed ASV table looks like this:

.. code-block:: bash

   run    marker    variant    sequence_length    read_count    tpos1_run1    tnegtag_run1    14ben01    14ben02    keep_tpos1_run1    clusterid    clustersize    chimera_borderline    sequence
   run1    MFZR    25    181    478    478    0    0    0    0    25    1    False    ACTATACCTTATCTTCGCAGTATTCTCAGGAATGCTAGGAACTGCTTTTAGTGTTCTTATTCGAATGGAACTAACATCTCCAGGTGTACAATACCTACAGGGAAACCACCAACTTTACAATGTAATCATTACAGCTCACGCATTCCTAATGATCTTTTTCATGGTTATGCCAGGACTTGTT
   run1    MFZR    51    181    165    0    0    0    165    0    51    1    False    ACTATATTTAATTTTTGCTGCAATTTCTGGTGTAGCAGGAACTACGCTTTCATTGTTTATTAGAGCTACATTAGCGACACCAAATTCTGGTGTTTTAGATTATAATTACCATTTGTATAATGTTATAGTTACGGGTCATGCTTTTTTGATGATCTTTTTTTTAGTAATGCCTGCTTTATTG
   run1    MFZR    88    175    640    640    0    0    0    1    88    1    False    ACTATATTTTATTTTTGGGGCTTGATCCGGAATGCTGGGCACCTCTCTAAGCCTTCTAATTCGTGCCGAGCTGGGGCACCCGGGTTCTTTAATTGGCGACGATCAAATTTACAATGTAATCGTCACAGCCCATGCTTTTATTATGATTTTTTTCATGGTTATGCCTATTATAATC
 
This step is particularly useful when analysing several runs with the same mock samples, since in the case the same *known_keep_occurrences.tsv* can be used for all runs.


.. _taxassign_reference:

The command taxassign
--------------------------------

All variants retained at the end of the filtering steps undergo a taxonomic assignation with this command. 

A quick introduction to the **vtam taxassign** command is given in the tutorial :ref:`taxassign_tutorial`. 

The command line arguments of the **taxassign** command can be obtained with 

.. code-block:: bash

    vtam taxassign --help

The most important arguments are these inputs and outputs:

Inputs:
    - :ref:`asvtable <ASVtable_io>`: TSV file containing variants in the last column (using ‘sequence’ as a header). Typically it is an ASV table. 
    - **sorteddir**: Directory containing the demultiplexed fasta files.
Database:
    - :ref:`db <db_io>`: Name of an sqlite db. It is used to store and access information on runs, markers, samples, variants, filtering steps and taxonomic assignations. It is created by the **filter** command if it does not exist and it is completed if it exists already. 
    - :ref:`taxonomy <taxonomy_io>`: TSV file containing taxonomic information <LINK to the taxonomy.tsv>
    - :ref:`blastdbdir <BLAST_database_reference>`: directory containing the BLAST database <LINK to BLAST batabase>
    - :ref:`blastdbname <BLAST_database_reference>`: name of the BLAST database <LINK to BLAST batabase>
Outputs:
    - :ref:`output <output_io>`: The input TSV file completed by taxonomic assignations

Variants are BLASTed against the NCBI nt or custom BLAST database. Taxa name and rank is chosen based on lineages of the best hits. Lineages are constructed based on a taxonomy TSV file (See :ref:`Reference section <BLAST_database_reference>`).

For a given %identity between the variant and the :ref:`hits <BLASThit_glossary>`, select hits with :ref:`coverage <coverage_glossary>` >= **min_query_coverage**. Depending on the %identity and the **ltg_rule_threshold**, there are two possible rules:

    #. %identity>= **ltg_rule_threshold**: Determine the Lowest Taxonomic Group (LTG)<link to the glossary>: Take the Lowest Taxonomic Group that contains at least **include_prop** % of the hits. Otherwise the LTG is not inferred at that %identity level.
    #. %identity< **ltg_rule_threshold**: Determine the Lowest Taxonomic Group (LTG)<link to the glossary> (using the same rule as previously) only if selected hits contain at least **min_number_of_taxa** different taxa. Otherwise the LTG is not inferred at that %identity level.

VTAM intends to establish LTG by using first a high %identity (100%). If LTG cannot be defined with this %identity, the %identity is decreased gradually (100%, 99%, 97%, 95%, 90%, 85%, 80%, 75%, 70%) till an LTG can be established. For each variant, the %identity used to infer the LTG is also given in the output. These values should not be ignored. If LTG could only be defined at 80% or lower identity level, the results are not very robust. 

Let’s see two examples using default values:

Example 1 of taxonomic assignation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A variant has produced 47 hits to *Beatis rhodani* sequences and 1 to a *Baetis cf. rhodani* sequence at 100 %identity. The ltg_rule_threshold is 97%. 

+----------------------+------------+------------+------------+------------+ 
|Taxon                 |Nb of hits  |%identity   |%coverage   |TaxID       |
+======================+============+============+============+============+ 
|*Baetis rhodani*      |47          |100         |   100      |189839      |
+----------------------+------------+------------+------------+------------+ 
|*Baetis cf. rhodani*  |1           |100         |   100      |1469487     |
+----------------------+------------+------------+------------+------------+ 

For all these alignments, the coverage was above 80 (**min_query_coverage**), so all of them are included in the assignment process.
LTG can be determined at the 100% identity level. Since 47 out of 48 sequences are from Baetis rhodani (47/48 > 90 (**include_prop**)), the LTG is *Baetis rhodani*.

Example 2 of taxonomic assignation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A variant has produced 11 hits to 4 taxa at 91-97 %identity. The ltg_rule_threshold is 97%.

+------------------------------+------------+----------+----------+-------+--------------------------------------------------------------+ 
|Taxon                         |Nb of hits  |%identity |%coverage |TaxID  |Lineage                                                       |
+==============================+============+==========+==========+=======+==============================================================+  
|*Procloeon pulchrum*          |2           |96        |   90     |1592912|...Ephemeroptera; Pisciforma; Baetidae; Procloeon             |
+------------------------------+------------+----------+----------+-------+--------------------------------------------------------------+ 
|*Procloeon rivulare*          |6           |90-93     |   100    |603535 |...Ephemeroptera; Pisciforma; Baetidae; Procloeon             |
+------------------------------+------------+----------+----------+-------+--------------------------------------------------------------+ 
|*Baetidae sp. BOLD:ACL6167*   |1           |91        |   100    |1808092|... Ephemeroptera; Pisciforma; Baetidae; unclassified Baetidae|
+------------------------------+------------+----------+----------+-------+--------------------------------------------------------------+ 
|*Procloeon sp. BOLD:AAK9569*  |2           |91        |   100    |1712756|...Ephemeroptera; Pisciforma; Baetidae; Procloeon             |
+------------------------------+------------+----------+----------+-------+--------------------------------------------------------------+ 

For this variant there are no hits with at least 100% or 97% identity. 

There are only 2 hits with at least 95% identity, both coming from the same taxon. Since 95% is below the **ltg_rule_threshold** (97%), at this identity level we need at least 3 different taxa (**min_number_of_taxa**) to derive an LTG. Therefore, VTAM will not infer LTG at the 95% identity level. The idea of not inferring LTG if the % of identity and the number of taxa is low is to avoid assignment of a variant to taxonomic level which is too precise, and probably erroneous.

At 90% identity level, there are 11 hits, coming from 4 taxa (>min_number_of_taxa). All hits have at least 80% query coverage (**min_query_coverage**). All conditions are met to infer LTG, which is Procloeon. 

The *Baetidae sp. BOLD:ACL6167* is not included in the LTG, since the other sequences that are all part of *Procloeon* genus represent 10/11>90% (**include_prop**) of the hits. The idea of using **include_prop** as a parameter is to try to avoid partially or erroneously annotated sequences if there is sufficient data to support an assignment with a higher resolution. In this example **Baetidae sp. BOLD:ACL6167** can in fact be a **Procloeon** species, but we do not have this information.

These parameters can be passed in a YML file via the **params** argument (See :ref:`Numerical parameter file <numerical_parameter_file_reference>` ).

.. _optimize_reference:

The command optimize
--------------------------------

The optimization step aims to find the optimal values for the parameters of the filtering steps. It is based on occurrences that are clearly erroneous (false positives, flagged as ‘delete’) or clearly expected (flagged as ‘keep’). These occurrences are determined by the user. The optimization step will suggest a parameter setting that keeps all occurrences in the dataset flagged as ‘keep’ and delete most occurrences flagged as ‘delete’. A list of known ‘keep’ and ‘delete’ occurrences is given in the :ref:`kown_occurrences.tsv <known_occurrences_io>` input file. The preparation of this file can be facilitated by running the filter command with the **known_occurrences** option (:ref:`See Make the ASV table section <MakeAsvTable_reference>`).

**The ‘keep’ occurrences are the expected variants in the mock samples**. 

Delete occurrences are the following :
    - All occurences in negative controls. 
    - Unexpected occurrences in mock samples.
    - Some occurrences of variants in real samples can also be flagged ‘delete’ if there are samples from markedly different habitats (e.g. marine vs. freshwater; the presence of a variant assigned to a freshwater taxon in a marine sample is clearly an error, and the occurrence can be flagged as ‘delete’).

Running first the **filter** command with default parameters produces an ASV table, where most of the sequence artefacts are eliminated, and therefore the ASV table has a reasonable size to deal with in a spreadsheet. Users should identify ‘keep’ and ‘delete’ occurrences from this output before the optimization step.

All optimization steps are run on the original, non-filtered read counts. 

A quick introduction to the **optimize** command is given in the tutorial :ref:`optimize_tutorial`.

The command line arguments of the **optimize** command can be obtained with 

.. code-block:: bash

    vtam optimize --help

These arguments are the most important inputs and outputs:

Inputs:
    - :ref:`sortedinfo <sortedinfo_io>`: TSV with the information about each sorted read file (output of “sortreads“). The content of the file will define the dataset to be used by “optimize“.
    - **sorteddir**: Directory containing the demultiplexed fasta files.
    - :ref:`known_occurrences <known_occurrences_io>`: User created TSV file with occurrences clearly identified as 'delete', 'keep or 'tolerate' (See manual).
Database:
    - :ref:`db <db_io>`: Name of an sqlite db. It is used to store and access information on runs, markers, samples, variants, filtering steps and taxonomic assignations.
Outputs:
    - **outdir**: Path to a directory that will contain the following files:

The following files will be written in the "--outdir" directory <link>
    - :ref:`optimize_lfn_sample_replicate.tsv <optimize_lfn_sample_replicate_io>`
    - :ref:`optimize_lfn_read_count_and_lfn_variant.tsv OR optimize_lfn_read_count_and_lfn_variant_replicate.tsv<optimize_lfn_read_count_and_lfn_variant_io>`
    - :ref:`optimize_lfn_variant_specific.tsv OR optimize_lfn_variant_replicate_specific.tsv<optimize_lfn_variant_specific_io>`
    - :ref:`optimize_pcr_error.tsv <optimize_pcr_error_io>`


.. _OptimizePCRError_reference:

The step OptimizePCRError
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Optimizing the :ref:`pcr_error_var_prop <FilterPCRerror_reference>` parameter is based only on the mock sample composition. The idea is to detect unexpected variants highly similar to expected ones within the same sample. The unexpected variant is likely to come from a PCR error and the proportion of their read counts help to choose the value of **pcr_error_var_prop**.

For each mock sample, all occurrences of unexpected variants (all but the ‘keep’) with a single difference to an expected variant (‘keep’) are detected. For each of these unexpected occurrences *N(i_unexpected)j*/*N(i expected)j* is calculated. The value for **pcr_error_var_prop** should be higher than the maximum of these proportions. The results are sorted by run, marker, and then by *N_ij_unexpected_to_expected_ratio* in decreasing order. Therefore, for each run marker combination, the **pcr_error_var_prop** parameter should be set above the first value.

**Example of** *optimize_pcr_error.tsv* :

.. code-block:: bash

    run    marker    sample    variant_expected    N_ij_expected    variant_unexpected    N_ij_unexpected    N_ij_unexpected_to_expected_ratio    sequence_expected    sequence_unexpected
    run1    MFZR    tpos1_run1    264    3471    1051    63    0.01815039    ACTTTATTTTATTTTTGGTGCTTGATCAGGAATAGTAGGAACTTCTTTAAGAATTCTAATTCGAGCTGAATTAGGTCATGCCGGTTCATTAATTGGAGATGATCAAATTTATAATGTAATTGTAACTGCTCATGCTTTTGTAATAATTTTCTTTATAGTTATACCTATTTTAATT    CCTTTATTTTATTTTTGGTGCTTGATCAGGAATAGTAGGAACTTCTTTAAGAATTCTAATTCGAGCTGAATTAGGTCATGCCGGTTCATTAATTGGAGATGATCAAATTTATAATGTAATTGTAACTGCTCATGCTTTTGTAATAATTTTCTTTATAGTTATACCTATTTTAATT
    run1    MFZR    tpos1_run1    88    640    89    8    0.01250000    ACTATATTTTATTTTTGGGGCTTGATCCGGAATGCTGGGCACCTCTCTAAGCCTTCTAATTCGTGCCGAGCTGGGGCACCCGGGTTCTTTAATTGGCGACGATCAAATTTACAATGTAATCGTCACAGCCCATGCTTTTATTATGATTTTTTTCATGGTTATGCCTATTATAATC    ACTATATTTTATTTTTGGGGCTTGATCCGGAATGCTGGGCACCTCTCTAAGCCTTCTAATTCGTGCCGAGCTGGGGCACCCGGGTTCTTTAATTGGCGACGATCAAATTTACAATGTAATCGTCACAGCCCATGCTTTTATTATGATTTTTTTCATGGTTATGCCTATTATAATT


.. _OptimizeLFNsampleReplicate_reference:

The step OptimizeLFNsampleReplicate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Optimizing the :ref:`lfn_sample_replicate_cutoff <FilterLFN_reference>` parameter is based only on the mock sample composition.

For each ‘keep’ variant in all mock samples (*j*), in all replicates (*k*), *N_ijk*/*N_jk* is calculated and these proportions are ordered increasingly. (*N_jk* is the total number of reads in the sample *j* replikate *k*. It includes all variants, not just the ones marked as keep) The optimal value for **lfn_sample_replicate_cutoff**, should be under the smallest proportion. This ensures keeping of all expected variants by the **LFN_sample_replicate** filter. 

In a fairly complex mock sample this smallest proportion rarely exceeds 0.01, but it is more often in the order of 0.001. Avoid setting this value higher than 0.01, since it can eliminate too many real occurrences.

One scenario (quite rare in our experience) is that all keep variants amplify well, and the suggested value for **lfn_sample_replicate_cutoff** is relatively high (>0.01). If the real samples are expected to have a high number of taxa (e.g. > 20), it is better to lower this value (default in 0.001).

Another case is that one or several expected variants have very few reads compared to the total number of reads in a sample-replicate (*N_ijk*/*N_jk*), while for the majority of the replicates of the same sample *N_ijk*/*N_jk* is not too low. In the example below, variant 75 in sample Tpos2_prerun and replicate 2 has only 0.01% of the reads, while in the other 2 replicates it is higher than 0.4%. In this case it is better to choose 0.004 as a **lfn_sample_replicate_cutoff**, knowing that the variant is not filtered out in 2 out of the 3 replicates, therefore, it will remain in the sample after pooling the results over replicates. 
(See :ref:`FilterMinReplicateNumber <FilterMinReplicateNumber_reference>`)

**Example of** *optimize_lfn_sample_replicate.tsv* :

.. code-block:: bash

    run    marker    sample    replicate    variant_id    N_ijk    N_jk    lfn_sample_replicate: N_ijk/N_jk    round_down    sequence
    prerun    MFZR    Tpos2_prerun    2    75    1    10234    0.00010234    0.000100000    ACTATATTTTATTTTTGGGGCTTGATCCGGAATGCTGGGCACCTCTCTAAGCCTTCTAATTCGTGCCGAGCTGGGGCACCCGGGTTCTTTAATTGGCGACGATCAAATTTACAATGTAATCGTCACAGCCCATGCTTTTATTATGATTTTTTTCATGGTTATGCCTATTATAATC
    prerun    MFZR    Tpos1_prerun    3    75    42    10414    0.00403303    0.00400000    ACTATATTTTATTTTTGGGGCTTGATCCGGAATGCTGGGCACCTCTCTAAGCCTTCTAATTCGTGCCGAGCTGGGGCACCCGGGTTCTTTAATTGGCGACGATCAAATTTACAATGTAATCGTCACAGCCCATGCTTTTATTATGATTTTTTTCATGGTTATGCCTATTATAATC
    prerun    MFZR    Tpos2_prerun    3    75    63    15038    0.00418939    0.00410000


.. _OptimizeLFNReadCountAndLFNvariant_reference:

The steps OptimizeLFNReadCountAndLFNvariant and OptimizeLFNReadCountAndLFNvariantReplicate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These steps find the combinations **lfn_variant_cutoff/read_count_cutoff** and **lfn_variant_replicate_cutoff/read_count_cutoff** that minimize the number of ‘delete’ occurrences while keeping all ‘keep’ occurrences.

The **fitler_lfn_variant** and **filter_lfn_variant_replicate** are alternatives around the same idea: Filtering occurrences in function of the their read count in the sample-replicate compared to the total number of reads of the variant in the run (*N_ijk*/*N_i*; **filter_lfn_variant**) or in the replicate (*N_ijk*/*N_ik*; **filter_lfn_variant_replicate**). The command **optimize** can have **lfn_variant** or **lfn_variant replicate** mode. Just like for the **filter** command, the default is the **lfn_variant** mode and the **lfn_variant_replicate** mode can be activated by the **lfn_variant_replicate** flag in the command line (see <link>). For simplicity, we will use **filter_lfn_variant** in the rest of this section.

All **FilterLFN** steps and **FilterMinReplicateNumber** are run on the original non-filtered data using a large number of combinations of **lfn_variant_cutoff** and **read_count_cutoff** (all other parameters are default). The values for these two thresholds vary between their default value till the highest value that keeps all ‘keep’ occurrences. For each combination, the number of ‘delete’ occurrences remaining in the dataset are counted (nb_delete) and printed to a spreadsheet in increasing order. Users should choose the parameter combination with lowest nb_delete.

**Example of** *optimize_lfn_read_count_and_lfn_variant.tsv*:

.. code-block:: bash

    occurrence_nb_keep    occurrence_nb_delete    lfn_nijk_cutoff    lfn_variant_cutoff    run    marker
    6    4    73    0.001    run1    MFZR
    6    4    73    0.005    run1    MFZR


.. _pool_reference:

The command pool
--------------------------------

This command will pool the results of several run-marker combinations into one ASV table. Variants identical on their identical regions are regrouped to the same line.

Inputs:
    - :ref:`run_marker <runmarker_io>`: TSV file listing all run marker combinations to be pooled
Database:
    - :ref:`db <db_io>`: Name of an sqlite db. It is used to store and access information on runs, markers, samples, variants, filtering steps and taxonomic assignations. It is created by **filter** if it does not exist and it is completed if it exists already. 
Outputs:
    - :ref:`asvtable <asvtable_io>`: Name of the pooled ASV table

In order to increase the chance of amplifying a large number of taxa, more than one primer pair (markers) can be used to amplify the same locus. The annealing sites of different markers can be slightly different, making the direct pooling of the results of different markers impossible. This step is carried out using the **vsearch --cluster_size** command with 1 as identity cutoff. In this way only variants identical in their overlapping regions are pooled together. 

As for making ASV tables after filtering, all variants are also clustered using a user defined identity cutoff (**cluster_identity**; 097 by default). The **clusterid** is the centroid, and the number of variants are given in the **clustersize** column. This clustering is a simple help for the users to identify similar variants, but it is NOT part of the filtering process nor part of pooling identical variants on their overlapping regions. Clusters of different ASV tables are not directly comparable. 

The output table contains the following columns:

    - variant_id: ID of a variant representative of the variants identical on their overlapping region
    - pooled_variants: List of variant identical on their overlapping region
    - run
    - marker
    - [one column per sample] : presence(1)/absence(0)
    - clusterid: Certoïd of a cluster (cluster_identity threshold)
    - clustersize: Number of variants in the cluster
    - pooled_sequences: Comma separated list of variants identical on their overlapping region
    - sequence: Sequence of the representative variant


.. _taxonomy_reference:

The command taxonomy and the taxonomic lineage input
----------------------------------------------------------------

VTAM requires a taxonomical file in TSV format that looks like this:

.. code-block:: bash

    tax_id    parent_tax_id    rank    name_txt    old_tax_id
    1    1    no rank    root    
    2    131567    superkingdom    Bacteria    
    6    335928    genus    Azorhizobium    
    7    6    species    Azorhizobium caulinodans    395
    9    32199    species    Buchnera aphidicola    28241
    10    1706371    genus    Cellvibrio

These are the columns of the file:
    - tax_id: Integer with the NCBI taxonomic ID.
    - parent_tax_id: Integer with the parent’s NCBI taxonomic ID.
    - rank: String with the taxonomic rank.
    - name_txt: String with the scientific name.
    - old_tax_id: Optional. Integer with a previous taxonomic ID that will be tried if an entry was not found in the tax_id column.

A precomputed file can be downloaded using this command. This version in not necessarily the most up to date compared to the ncbi taxonomy database, but it works with our custom database:

.. code-block:: bash

    vtam taxonomy --output taxonomy.tsv --precomputed

Alternatively, the command will download the up-to-date ncbi taxonomy database (`<https://www.ncbi.nlm.nih.gov/taxonomy>`_) and create a fresh TSV file with the latest data in NCBI.
This is strongly recommended if you are using a recently downloaded version of the ncbi_nt:

.. code-block:: bash

    vtam taxonomy --output taxonomy.tsv

This step can take several minutes. Make sure you have a steady internet connection. 

The *taxonomy.tsv* file created by the script is ready to use as is if you intend to use the full NCBI nucleotide :ref:`BLAST database <BLAST_database_reference>` or our :ref:`precomputed non-redundant database specific to COI <non-redundant_COI_reference>`.  However, if you create a custom database, containing sequences from taxa not yet included in NCBI taxonomic database, you have to complete this file with arbitrary Taxonomic IDs used in your custom database and link them to existing NCBI taxids. We suggest using negative TaxIDs for taxa not present in NCBI Taxonomy database.

An example is found here:

.. code-block:: bash

    tax_id    parent_tax_id    rank    name_txt    old_tax_id
    ...
    233297    34643    genus    Chrysso    449633
    41073    535382    family    Carabidae 
    ...
    -1    233297    species    Chrysso pelyx
    -2    41073    genus    Pelophila
    -3    -2    species    Pelophila borealis


.. _BLAST_database_reference:

The BLAST database
----------------------------------------------------------------

VTAM uses a BLAST database for taxonomic assignment. There are three possibilities:

    - Download the full NCBI nucleotide (NCBI_nt) BLAST (ca. 100 Gb).
    - Download a precomputed non-redundant database specific to the first half (600-700bp) of the COI gene (less than 1 Gb).
    - Build a custom BLAST database with local sequences.

.. _NCBI_nt_reference:

NCBI nt BLAST database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The NCBI nucleotide (Genbank) database can be downloaded with the following commands:

.. code-block:: bash

    wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt*"
    md5sum --check nt.*.tar.gz.md5 
    tar -zxvf nt.*.tar.gz 
    rm nt.*.tar.gz*

It contains all NCBI nucleotide sequences and thus, it is not limited to COI.

.. _non-redundant_COI_reference:

Precomputed non-redundant COI database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We have created a non-redundant database specific to the first half (600-700bp) of the COI gene. It contains COI sequences from NCBI nt and BOLD, they cover at least 80% of the barcoding region of the COI. Identical sequences of the same taxon are present only once in the database. The latest version of this database can be downloaded with the following command:

.. code-block:: bash

    vtam coi_blast_db --blastdbdir vtam_db/coi_blast_db

Earlier versions can be downloaded by specifying the name of database:

.. code-block:: bash

    vtam coi_blast_db --blastdbdir vtam_db/coi_blast_db --blastdbname coi_blast_db_20200420

The available versions are found here: `<https://github.com/aitgon/vtam/releases/latest>`_


Custom database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To create a custom database, make sure that you have installed blast>=2.9.0 available in Conda.

Then create a fasta file *mydb.fas* and a mapping *taxid_map.tsv* from sequences to NCBI taxa. These are the first lines of the *mydb.fas* file.

.. code-block:: bash

	>MG1234511
	AACACTTTATTTTATTTTTGGAATTTGAGCTGGAATAGTAGGAACATCATTAAGAATTTTAATTCGATTGGAATTAAGAACAATTTCTAATTTAATTGGAAATGATCAAAATGTAATTGTAACAGNNNCTTTTATTATAATTTTCTTTATAGTTATACCTATTTTAATT
	>MG1588471_1
	AACATTATATTTTATTTTTGGTGCTTGATCAAGAATAGTGGGGACTTCTTTAAGAATACTTATTCGAGCTGAATTAGGGTGTCCGAATGCTTTAATTGGGGATGACCAAATTTATAATGTTATTGTAACTGCCCATGCTTTTATTATAAttttttttATAGTAATACCTATTATAATT
	>MG1588471_2
	AACATTATATTTTATTTTTGGTGCTTGATCAAGAATAGTGGGGACTTCTTTAAGAATACTTATTCGAGCTGAATTAGGGTGTCCGAATGCTTTAATTGGGGATGACCAAATTTATAATGTTATTGTAACTGCCCATGCTTTTATTATAAttttttttATAGTAATACCTATTATAATT
	>LOCAL_SEQ1
	TATTTTATTTTTGGAATATGAGCAGGAATATTAGGATCATCAATAAGATTAATTATTCGAATAGAACTAGGTAACCCTGGATTTTTAATTAATAATGATCAAATTTACAATTCTTTTGTAACAGCTCATGCATTTATTATAAttttttttATAGTAATACCAATTATAATT

These are the first lines of the *taxid_map.tsv* file where the first column is the sequence Identifier and the second is the NCBI TaxID of the taxon of origin or arbitrary unique TaxIDs included in the :ref:`taxonomy.tsv  <taxonomy_reference>` file:

.. code-block:: bash

	MG1234511    2384799
	MG1588471_1    2416844
	MG1588471_2    2416875
	LOCAL_SEQ1    -3

Then you run this command:

.. code-block:: bash

    makeblastdb -in mydb.fas -parse_seqids -dbtype nucl -out nt -taxid_map taxid_map.tsv


Traceability
----------------------------------------------------------------

Variants, samples, the results of filtering steps, and the taxonomic assignments are stored in a sqlite database. This provides a possibility to trace the history of the analyses. You can easily discover the sqlite database with a sqlite browser program ( For example `<https://sqlitebrowser.org/>`_ or `<https://sqlitestudio.pl>`_).

Here we propose a few examples of sql commands that can be adapted to extract the information you need.

Each filter has a table in the database, and they contain the following fields:
    - run_id
    - marker_id
    - varinat_id
    - sample_id
    - replicate
    - read_count
    - filter_delete
    
The FilterLNF table is special, because it is composed of several filters: filter_id=2, 3, ... To select variants that passed all filters, we need to check filter_id=8:

Count the number of variants after FilterChimera for run 1 and marker 1

.. code-block:: bash

    select count(distinct variant_id) from FilterChimera where run_id=1 and marker_id=1 and filter_delete=0
    # for FilterLFN
    select count(distinct variant_id) from FilterLFN where run_id=1 and marker_id=1 and filter_delete=0 and filter_id=8

Count the number of non_empty samples after FilterChimera for run 1

.. code-block:: bash

    select count(distinct sample_id) from FilterChimera where run_id=1 and filter_delete=0
    # for FilterLFN
    select count(distinct sample_id) from FilterLFN where run_id=1 and filter_delete=0 and filter_id=8

Count the number of reads that passed the filter

.. code-block:: bash

    select sum(read_count) from FilterChimera where run_id=1 and filter_delete=0
    # for FilterLFN
    select sum(read_count) from FilterLFN where run_id=1 and filter_delete=0 and filter_id=8

Get the list of the samples after a given filtering step

.. code-block:: bash

    select distinct Sample.name as sample from FilterChimera, Sample where FilterChimera.sample_id=Sample.id and filter_delete=0 order by Sample.name
    # for FilterLFN
    select distinct Sample.name as sample from FilterLFN, Sample where FilterLFN.sample_id=Sample.id and filter_delete=0 and filter_id=8 order by Sample.name

Get the list of the sample-replicates after a given filtering step

.. code-block:: bash

    select distinct Sample.name as sample, FilterChimera.replicate from FilterChimera, Sample where FilterChimera.sample_id=Sample.id and filter_delete=0 order by Sample.name, FilterChimera.replicate
    # for FilterLFN
    select distinct Sample.name as sample, FilterLFN.replicate from FilterLFN, Sample where FilterLFN.sample_id=Sample.id and filter_delete=0 and filter_id=8 order by Sample.name, FilterLFN.replicate

Get the list of the variants after a given filtering step

.. code-block:: bash

    select distinct Variant.id, Variant.sequence from FilterChimera, Variant where FilterChimera.variant_id=Variant.id and filter_delete=0 order by Variant.id
    # for FilterLFN
    select distinct Variant.id, Variant.sequence from FilterLFN, Variant where FilterLFN.variant_id=Variant.id and filter_delete=0 and filter_id=8 order by Variant.id
