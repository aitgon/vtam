Developper
=================================================

Taxon assignation method
------------------------------------------------------

Inputs:
------------------------------------------------------
variant_sequence_fasta: Fasta file with 1 variant at a time.

database_fasta: Fasta file with all reference sequence to be aligned on the 
variant sequence.

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




	

	
