Developer
=================================================

Taxon assignation method
------------------------------------------------------

Input
------------------------------------------------------

- variant_sequence_fasta: Fasta file with 1 variant at a time.
- database_fasta: Fasta file with all reference sequence to be aligned on the variant sequence.
- tax assign parameters: TSV table file  with four columns

    * identity_threshold
    * min_tax_level
    * max_tax_resolution
    * min_tax_n

Example of tax assign parameters:

*tax_assign_pars.tsv*

.. code-block:: bash

    100.0	species	subspecies	1
    97.0	genus	species	1
    95.0	family	species	3
    90.0	order	family	3
    85.0	order	order	3
    80.0	class	order	5

Step 0: Database creation with database_fasta:
------------------------------------------------------
The sequences ids of database_fasta contains needed information as taxon id,
taxonomic rank etc ... This information will be stored in a database table for
be used on all variants and gain time.

Step1: Vsearch
------------------------------------------------------

Goal: Align the variant sequence on the database

.. code-block:: bash

    variant_sequence_fasta ---------
                                    \
                                     ---> vsearch --------> alignement_result.tsv
                                    /
    database_fasta -----------------

The output tsv must be on the following shape:

.. code-block:: bash

    TCTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGGACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT	6764813	100.0

Step 2: Taxonomic association parameters
------------------------------------------------------

The user have to give a tsv file with the following information:

.. code-block:: bash

    idx		min_taxon_level	max_taxon_resolution	min_taxon_n
    100.0	species	subspecies	1
    97.0	genus	species	1

idx: Is the minimum percentage "quality" of alignment required.
min_taxon_level: less detailed taxonomy level possible for idx percentage
max_taxon_resolution: more accurate taxonomy level accorded to a idx
min_taxon_n: minimum of hits after filter required for an taxonomic association

Step3: Looping on idx and 1st filters application
------------------------------------------------------

For each idx containend in taxassign parameters (Allow to use the join information as min_taxon_level, ...):

    - Selecting hits from the vsearch alignment (dataframe), with the idx. If the selected_dataframe is empty, pass to the next
      idx.

    - Remove hits from selected_dataframe when rank_name >= min_taxon_level. If the selected_dataframe is empty, pass to the next

    - If there are less lines remaining in the selected_dataframe. Pass to the next idx.

Step 4: Phylogenetic lineage creation
------------------------------------------------------

After removing non conform its.

We need to create the phylogenomic lineage for every hit.

For that we get all the taxon_seq_id (id of the hits), and we make a request on the database to get the following information:

.. code-block:: bash

    taxon_seq_id tax_name   tax_id  rank    parent_tax_id

    5244419 Echinorhynchida	57283   order	45080
    5244429	Echinorhynchida	57283	order	45080





Schema:

for each vsearch identity idx ()

    select vsearch_output_lines

    if idx(line) >= idx -------True------> line selected
                        \
                         ------False---- > line non selected | If selected_lines == 0: Continue with another idx

    remove from df_selected when rank_name >= min_taxon_level

    is_true_min_taxon_n -------False-----> next idx
                        \
                         ------True-------> Following instruction

                        for each df_selected_taxon_id
                            create phylogenetic_lineage

                        get lower taxonomic group (LTG)
                            with <indice_prop>*100% of the hits sequences

                        is_true_max_resolution ------True------> Keep LTG
                                               \
                                                ---------------> increase rank_name up to max_tax_resolution

At the end of this step, we get this data frame (*tax_lineage_df*) with these columns:

- tax_seq_id: Hits of vsearch aligned to variant with similarity percentage above a given threshold
- following columns 

.. code-block:: bash

    tax_seq_id	no rank	phylum	class	subclass	infraclass	superorder	order	suborder	infraorder	family	subfamily	genus	species
    6320345	131567	6073	6101	6102			6103	86626.0		37511.0		6115.0	86610
    4307609	131567	6073	6101	6102			6125						
    4314607	131567	6073	6101	6102			6125						
    2658650	131567	6656	50557	7496	33340.0		7041						
    2658649	131567	6656	50557	7496	33340.0		7041						
    6349460	131567	6073	6101	6102			6103	86626.0		37511.0		6115.0	86610
    6349457	131567	6073	6101	6102			6103	86626.0		37511.0		6115.0	86610
    8073839	131567	6656	50557	7496	33340.0		7147	7203.0	43733.0	7371.0	43914.0	7374.0	65466
    6297084	131567	6073	6101	6102			6103	86626.0		37511.0		6115.0	
    6349463	131567	6073	6101	6102			6103	86626.0		37511.0		6115.0	86610
    5144903	131567	6656	6854	6933		6934.0	34634	281668.0	1723665.0	99213.0			
    5285255	131567	6656	6854	6933		6934.0	34634						
    7492317	131567	6656	6854	6933		6934.0	34634	281668.0	1723665.0	99224.0		99225.0	
    6288281	131567	6073	6101	6102			6103	86626.0		37511.0		6115.0	
    6349397	131567	6073	6101	6102			6103	86626.0		37511.0		6115.0	6116

Step 5: LTG assignement
------------------------------------------------------

Given the taxon lineage data frame(*tax_lineage_df*), here we search for the low taxonomy group (LTG) with these rules

- The LTG must comprise a percentage (90%) of hits (Number of rows in *tax_lineage_df*)
- The rank of the LTG must be more detailed than *min_tax_level*

The more detailed taxon and its rank following these rule will be set as a temporary *ltg* and *ltg_rank*. Then two situations:

- *ltg_rank* less or equally detailed than *max_tax_resolution_id*. Then we keep this LTG
- *ltg_rank* more detailed than *max_tax_resolution*.
    * Then go up in the taxonomic line of current LTG up to *max_tax_resolution* rank of current LTG and set *max_tax_resolution* taxon to new LTG
        + If *max_tax_resolution* taxon is not defined, then keep current LTG as LTG

