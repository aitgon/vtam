import sqlite3

import pandas

from wopmetabarcoding.utils.constants import rank_hierarchy


def taxonomy_sqlite_to_df(taxonomy_sqlite):
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

def tax_id_to_taxonomy_lineage(tax_id, taxonomy_db_df):
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



def select_ltg(tax_lineage_df, identity, identity_threshold=97, include_prop=90, min_number_of_taxa=3):
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
        putative_ltg_df['putative_ltg_percentage'] = putative_ltg_df.putative_ltg_count / putative_ltg_df.shape[0] * 100
        ltg_tax_id = putative_ltg_df.loc[putative_ltg_df.putative_ltg_percentage >= include_prop, 'putative_ltg_id'].tail(1).values[0]
        ltg_rank = putative_ltg_df.loc[putative_ltg_df.putative_ltg_percentage >= include_prop, 'putative_ltg_id'].index[-1]
    else: # rule for min_number_of_taxa
        ltg_tax_id = int(putative_ltg_df.loc[putative_ltg_df.putative_ltg_count >= min_number_of_taxa, 'putative_ltg_id'].tail(1).values[0])
        ltg_rank = putative_ltg_df.loc[putative_ltg_df.putative_ltg_count >= min_number_of_taxa, 'putative_ltg_id'].index[-1]
    return ltg_tax_id, ltg_rank
