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

At the end of this step, we get this kind of data frame (*tax_lineage_df*).

.. code-block:: bash

    tax_seq_id	species	genus	subfamily	family	infraorder	suborder	order	superorder	infraclass	subclass	class	phylum	no rank
    6320345	86610	6115.0		37511.0		86626.0	6103			6102	6101	6073	131567
    4307609							6125			6102	6101	6073	131567
    4314607							6125			6102	6101	6073	131567
    2658650							7041		33340.0	7496	50557	6656	131567
    2658649							7041		33340.0	7496	50557	6656	131567
    6349460	86610	6115.0		37511.0		86626.0	6103			6102	6101	6073	131567
    6349457	86610	6115.0		37511.0		86626.0	6103			6102	6101	6073	131567
    8073839	65466	7374.0	43914.0	7371.0	43733.0	7203.0	7147		33340.0	7496	50557	6656	131567
    6297084		6115.0		37511.0		86626.0	6103			6102	6101	6073	131567
    6349463	86610	6115.0		37511.0		86626.0	6103			6102	6101	6073	131567
    5144903				99213.0	1723665.0	281668.0	34634	6934.0		6933	6854	6656	131567
    5285255							34634	6934.0		6933	6854	6656	131567
    7492317		99225.0		99224.0	1723665.0	281668.0	34634	6934.0		6933	6854	6656	131567
    6288281		6115.0		37511.0		86626.0	6103			6102	6101	6073	131567
    6349397	6116	6115.0		37511.0		86626.0	6103			6102	6101	6073	131567

Step 5: LTG assignement
------------------------------------------------------

Given the tax_lineage (*tax_lineage_df*), the low taxonomy group (LTG) is a tax_id, that follows these contraints

- must comprise a percentage (90%) of hits (Number of rows in *tax_lineage_df*)
- its rank must be more detailed than *min_tax_level*
- its rank must be less detailed than *max_tax_resolution*

