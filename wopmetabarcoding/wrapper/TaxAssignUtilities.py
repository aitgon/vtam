import os
import pandas
import sqlite3

from wopmetabarcoding.utils.constants import rank_hierarchy
from wopmetabarcoding.utils.logger import logger


def f01_taxonomy_sqlite_to_df(taxonomy_sqlite):
    """
    Imports taxonomy_db.sqlite file with taxonomy table (tax_id, parent_tax_id, rank, name_txt, old_tax_id)
    into DataFrame

    Args:
        taxonomy_sqlite (String): Path to taxonomy_db.sqlite


    Returns:
        pandas.DataFrame: with taxonomy information

    """
    con = sqlite3.connect(taxonomy_sqlite)
    taxonomy_db_df = pandas.read_sql_query("SELECT * FROM taxonomy", con=con)
    con.close()
    return taxonomy_db_df


def f02_variant_df_to_fasta(variant_df, fasta_path):
    """
    Takes variant DF with two columns (variant_id, variant_sequence) and return FASTA file path

    Args:
        variant_df (pandas.DataFrame): DF with two columns (variant_id, variant_sequence)
        fasta_path (str): Path to FASTA file


    Returns:
        None

    """
    with open(fasta_path, "w") as fout:
        for row in variant_df.itertuples():
            fout.write(">{}\n{}\n".format(row.variant_id, row.variant_sequence))


def f03_blast(variant_id, variant_sequence):
    """
    Runs Blast

    Args:
        variant_id (Integer): Internal ID of variant
        variant_sequence (String): Sequence of variant

    Returns:
        String: Path to output of blast in TSV format
    """
    # http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc98
    # Create FASTA
    from Bio.Blast.Applications import NcbiblastnCommandline
    blastn_cline = NcbiblastnCommandline(query="opuntia.fasta", db="nr", evalue=0.001, outfmt=5, out="opuntia.xml")


def f04_import_blast_output_into_df(blast_output_tsv_path):
    """
    Imports Blast output TSV into DataFrame

    Args:
        blast_output_tsv_path (String): Path to output of blast in TSV format

    Returns:
        DataFrame: with three columns: target_id, identity, target_tax_id
    """
    blast_out_df = pandas.read_csv(blast_output_tsv_path, sep="\t", header=None, names=['target_id', 'identity', 'target_tax_id'], usecols=[1, 2, 5])
    # extract the target_id
    blast_out_df.target_id = blast_out_df.target_id.str.split('|', 2).str[1]
    blast_out_df.target_id = pandas.to_numeric(blast_out_df.target_id)
    return blast_out_df


def f04_1_tax_id_to_taxonomy_lineage(tax_id, taxonomy_db_df):
    """
    Takes tax_id and taxonomy_db.sqlite DB and create a dictionary with the taxonomy lineage

    Args:
        tax_id (Integer): Identifier of taxon
        taxonomy_db_df (pandas.DataFrame: DataFrame with taxonomy information


    Returns:
        Dictionnary: with taxonomy information for given tax_id, {'tax_id': 183142, 'species': 183142, 'genus': 10194, 'family': 10193, 'order': 84394,
                         'superorder': 1709201, 'class': 10191, 'phylum': 10190, 'no rank': 131567, 'kingdom': 33208,
                         'superkingdom': 2759}

    """
    lineage_dic = {}
    lineage_dic['tax_id'] = tax_id
    while tax_id != 1:
        tax_id_row = taxonomy_db_df.loc[taxonomy_db_df.tax_id == tax_id,]
        rank = tax_id_row['rank'].values[0]
        parent_tax_id = tax_id_row['parent_tax_id'].values[0]
        lineage_dic[rank] = tax_id
        tax_id = parent_tax_id
    return lineage_dic

def f05_blast_result_subset(blast_result_subset_df, taxonomy_db_df):
    """
    Takes blast_result_subset_df and returns tax_lineage_df with lineages of each target_tax_id

    Args:
        blast_result_subset_df (DataFrame): DataFrame with result subset for given identity
            where each column is target_id and target_tax_id
    target_id  target_tax_id
0  1049499563         761875
1  1049496963         761875
        taxonomy_db_df (pandas.DataFrame: DataFrame with taxonomy information

    Returns:
        DataFrame: with a lineage per row that corresponds to the the lineage of each taxon from the blast sequence targets
   class  cohort  family  genus  infraclass  kingdom  no rank  order  phylum  species  subclass  subfamily  suborder  subphylum  superfamily  superkingdom  superorder
0  50557   33392   41030  50443       33340    33208   131567  30263    6656   761875      7496     147297     93873       6960        41029          2759       85604
1  50557   33392   41030  50443       33340    33208   131567  30263    6656   761875      7496     147297     93873       6960        41029          2759       85604

    """
    lineage_list = []
    for target_tax_id in blast_result_subset_df.target_tax_id.unique().tolist():
        lineage_list.append(f04_1_tax_id_to_taxonomy_lineage(target_tax_id, taxonomy_db_df))
    tax_lineage_df = pandas.DataFrame(lineage_list)
    tax_lineage_df = blast_result_subset_df.merge(tax_lineage_df, left_on='target_tax_id', right_on='tax_id')
    tax_lineage_df.drop('target_id', axis=1, inplace=True)
    tax_lineage_df.drop('target_tax_id', axis=1, inplace=True)
    return tax_lineage_df

def f06_select_ltg(tax_lineage_df, identity, identity_threshold=97, include_prop=90, min_number_of_taxa=3):
    """
    Given tax_lineage_df, selects the LTG

    Args:
        tax_lineage_df (pandas.DataFrame): DF where each column is a rank, rows are different target_ids and values are putative_ltg_ids
        identity (int): Identity value that were used to select the results from Blast
        identity_threshold (int): Identity value where we change of using include_prop method to min_number_of_taxa, default 97
        include_prop (int): Percentage out of total selected blast hits for LTG to be present when identity>=identity_threshold
        min_number_of_taxa (int): Minimal number of taxa, where LTF must be present when identity<identity_threshold

    Returns:
        ltg_tax_id (int): Taxonomical ID of LTG
        ltg_rank (str): Rank of LTG

    """
    lineage_list_df_columns_sorted = list(filter(lambda x: x in tax_lineage_df.columns.tolist(), rank_hierarchy))
    tax_lineage_df = tax_lineage_df[lineage_list_df_columns_sorted]
    putative_ltg_df = pandas.DataFrame(
        {'putative_ltg_id': tax_lineage_df.apply(lambda x: x.value_counts().index[0], axis=0),
         'putative_ltg_count': tax_lineage_df.apply(lambda x: x.value_counts().iloc[0], axis=0)})
    if identity >= identity_threshold: # rule for include_prop
        putative_ltg_df['putative_ltg_percentage'] = putative_ltg_df.putative_ltg_count / tax_lineage_df.shape[0] * 100
        ltg_tax_id = putative_ltg_df.loc[putative_ltg_df.putative_ltg_percentage >= include_prop, 'putative_ltg_id'].tail(1).values[0]
        ltg_rank = putative_ltg_df.loc[putative_ltg_df.putative_ltg_percentage >= include_prop, 'putative_ltg_id'].index[-1]
    else: # rule for min_number_of_taxa
        ltg_tax_id = int(putative_ltg_df.loc[putative_ltg_df.putative_ltg_count >= min_number_of_taxa, 'putative_ltg_id'].tail(1).values[0])
        ltg_rank = putative_ltg_df.loc[putative_ltg_df.putative_ltg_count >= min_number_of_taxa, 'putative_ltg_id'].index[-1]
    return ltg_tax_id, ltg_rank
