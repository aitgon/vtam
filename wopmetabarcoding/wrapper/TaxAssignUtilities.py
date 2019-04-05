import sqlite3

import pandas


def taxonomy_sqlite_to_df(taxonomy_sqlite):
    con = sqlite3.connect(taxonomy_sqlite)
    taxonomy_db_df = pandas.read_sql_query("SELECT * FROM taxonomy", con=con)
    con.close()
    return taxonomy_db_df

def tax_id_to_taxonomy_lineage(tax_id, taxonomy_db_df):
    lineage_dic = {}
    lineage_dic['tax_id'] = tax_id
    while tax_id != 1:
        tax_id_row = taxonomy_db_df.loc[taxonomy_db_df.tax_id == tax_id,]
        rank = tax_id_row['rank'].values[0]
        parent_tax_id = tax_id_row['parent_tax_id'].values[0]
        lineage_dic[rank] = tax_id
        tax_id = parent_tax_id
    return lineage_dic

