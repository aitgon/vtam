import os
import pandas
import sqlite3

from vtam.utils.constants import rank_hierarchy

from vtam.utils.constants import identity_list

from Bio.Blast.Applications import NcbiblastnCommandline

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
    Takes variant DF with two columns (variant_id, variant_sequence) and return FASTA file output

    Args:
        variant_df (pandas.DataFrame): DF with two columns (variant_id, variant_sequence)
        fasta_path (str): Path to FASTA file


    Returns:
        None

    """
    with open(fasta_path, "w") as fout:
        for row in variant_df.itertuples():
            fout.write(">{}\n{}\n".format(row.variant_id, row.variant_sequence))


def f03_lblast(variant_fasta_path, variant_tsv_path):
    """
    Runs Blast

    Args:
        variant_id (Integer): Internal ID of variant
        variant_sequence (String): Sequence of variant

    Returns:
        String: Path to output of qblast in TSV format
    """
    # http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc98
    # Create FASTA
    blastn_cline = NcbiblastnCommandline(query=variant_fasta_path, db="nt", evalue=1e-5,
                                         outfmt='"6 qseqid sacc pident evalue qcovhsp staxids"', dust='yes',
                                         qcov_hsp_perc=80, num_threads=1, out=variant_tsv_path)
    stdout, stderr = blastn_cline()


def f04_import_qblast_output_into_df(qblast_output_tsv_path):
    """
    Imports Blast output TSV into DataFrame

    Args:
        qblast_output_tsv_path (String): Path to output of qblast in TSV format

    Returns:
        DataFrame: with three columns: target_id, identity, target_tax_id
    """
    qblast_out_df = pandas.read_csv(qblast_output_tsv_path, sep="\t", header=None, names=['target_id', 'identity', 'target_tax_id'], usecols=[1, 2, 5])
    # extract the target_id
    qblast_out_df.target_id = qblast_out_df.target_id.str.split('|', 2).str[1]
    qblast_out_df.target_id = pandas.to_numeric(qblast_out_df.target_id)
    return qblast_out_df


def f04_1_tax_id_to_taxonomy_lineage(tax_id, taxonomy_db_df, give_tax_name=False):
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
        # try to use taxonomy_db_df.tax_id
        tax_id_row = taxonomy_db_df.loc[taxonomy_db_df.tax_id == tax_id,]
        # row empty, try to use old_tax_id
        if tax_id_row.shape[0] == 0:
            tax_id_row = taxonomy_db_df.loc[taxonomy_db_df.old_tax_id == tax_id,]
        rank = tax_id_row['rank'].values[0]
        parent_tax_id = tax_id_row['parent_tax_id'].values[0]
        lineage_dic[rank] = tax_id
        if give_tax_name: # return tax_name instead of tax_d
            tax_name = tax_id_row['name_txt'].values[0]
            lineage_dic[rank] = tax_name
        tax_id = parent_tax_id
    return lineage_dic

def f05_blast_result_subset(qblast_result_subset_df, taxonomy_db_df):
    """
    Takes qblast_result_subset_df and returns tax_lineage_df with lineages of each target_tax_id

    Args:
        qblast_result_subset_df (DataFrame): DataFrame with result subset for given identity
            where each column is target_id and target_tax_id
    target_id  target_tax_id
0  1049499563         761875
1  1049496963         761875
        taxonomy_db_df (pandas.DataFrame: DataFrame with taxonomy information

    Returns:
        DataFrame: with a lineage per row that corresponds to the the lineage of each taxon from the qblast sequence targets
   class  cohort  family  genus  infraclass  kingdom  no rank  order  phylum  species  subclass  subfamily  suborder  subphylum  superfamily  superkingdom  superorder
0  50557   33392   41030  50443       33340    33208   131567  30263    6656   761875      7496     147297     93873       6960        41029          2759       85604
1  50557   33392   41030  50443       33340    33208   131567  30263    6656   761875      7496     147297     93873       6960        41029          2759       85604

    """
    lineage_list = []
    for target_tax_id in qblast_result_subset_df.target_tax_id.unique().tolist():
        lineage_list.append(f04_1_tax_id_to_taxonomy_lineage(target_tax_id, taxonomy_db_df))
    tax_lineage_df = pandas.DataFrame(lineage_list)
    tax_lineage_df = qblast_result_subset_df.merge(tax_lineage_df, left_on='target_tax_id', right_on='tax_id')
    tax_lineage_df.drop('target_id', axis=1, inplace=True)
    tax_lineage_df.drop('target_tax_id', axis=1, inplace=True)
    return tax_lineage_df

# def f06_select_ltg(tax_lineage_df, identity, identity_threshold, include_prop, min_number_of_taxa):
def f06_select_ltg(tax_lineage_df, include_prop):
    """
    Given tax_lineage_df, selects the Ltg

    Args:
        tax_lineage_df (pandas.DataFrame): DF where each column is a rank, rows are different target_ids and values are putative_ltg_ids
        include_prop (int): Percentage out of total selected qblast hits for Ltg to be present when identity>=identity_threshold

    Returns:
        ltg_tax_id (int): Taxonomical ID of Ltg
        ltg_rank (str): Rank of Ltg

    """
    # Remove column if all values=NaN
    tax_lineage_df.dropna(axis='columns', how='all', inplace=True)
    lineage_list_df_columns_sorted = list(filter(lambda x: x in tax_lineage_df.columns.tolist(), rank_hierarchy))
    tax_lineage_df = tax_lineage_df[lineage_list_df_columns_sorted]
    putative_ltg_df = pandas.DataFrame({'putative_ltg_id': tax_lineage_df.apply(lambda x: x.value_counts().index[0], axis=0),'putative_ltg_count': tax_lineage_df.apply(lambda x: x.value_counts().iloc[0], axis=0)})
    # if identity >= identity_threshold: # rule for include_prop
    putative_ltg_df['putative_ltg_percentage'] = putative_ltg_df.putative_ltg_count / tax_lineage_df.shape[0] * 100
    ltg_tax_id = putative_ltg_df.loc[putative_ltg_df.putative_ltg_percentage >= include_prop, 'putative_ltg_id'].tail(1).values[0]
    ltg_rank = putative_ltg_df.loc[putative_ltg_df.putative_ltg_percentage >= include_prop, 'putative_ltg_id'].index[-1]
    # else: # rule for min_number_of_taxa
    #     ltg_tax_id = int(putative_ltg_df.loc[putative_ltg_df.putative_ltg_count >= min_number_of_taxa, 'putative_ltg_id'].tail(1).values[0])
    #     ltg_rank = putative_ltg_df.loc[putative_ltg_df.putative_ltg_count >= min_number_of_taxa, 'putative_ltg_id'].index[-1]
    return ltg_tax_id, ltg_rank

def f07_blast_result_to_ltg_tax_id(variantid_identity_lineage_df, identity_threshold, include_prop, min_number_of_taxa):
    """
    Main function that takes blast result with variant_id, target_id, identity and tax_id and returns ltg_tax_id and ltg_rank

    Example of the input variantid_identity_lineage_df
    variant_id  identity  target_tax_id  class  family   genus  infraclass  kingdom  no rank  order  phylum  species  subclass  suborder  subphylum  superkingdom
0           3   100.000         189839  50557  172515  189838       33339    33208   131567  30073    6656   189839      7496    185809       6960          2759
1           3    99.429         189839  50557  172515  189838       33339    33208   131567  30073    6656   189839      7496    185809       6960          2759
2           3    99.429         189839  50557  172515  189838       33339    33208   131567  30073    6656   189839      7496    185809       6960          2759
3           3    99.429         189839  50557  172515  189838       33339    33208   131567  30073    6656   189839      7496    185809       6960          2759
4           3    99.429         189839  50557  172515  189838       33339    33208   131567  30073    6656   189839      7496    185809       6960          2759
    Example of the output df:
    identity ltg_rank  ltg_tax_id  variant_id
0       100  species      189839           3
1       100  species     1077837           7
2        99  species     1077837           9

    Args:
        variantid_identity_lineage_df (pandas.DataFrame): DF with columns: variant_id, identity, target_tax_id and lineage_columns.
        identity_threshold (int): Identity value where we change of using include_prop method to min_number_of_taxa, default 97
        include_prop (int): Percentage out of total selected qblast hits for Ltg to be present when identity>=identity_threshold
        min_number_of_taxa (int): Minimal number of taxa, where LTF must be present when identity<identity_threshold

    Returns:
        ltg_df (pandas.DataFrame): DF with variant_id, ltg_tax_id and ltg_tax_rank

    """
    #
    list_variant_id_to_ltg = []
    variant_id_list = sorted(variantid_identity_lineage_df.variant_id.unique().tolist())
    #
    for variant_id in variant_id_list:  #  Loop sorted each variant
        for identity in identity_list:  # For each variant, loop each decreasing identity
            tax_lineage_by_variant_id_df = variantid_identity_lineage_df.loc[
                ((variantid_identity_lineage_df['variant_id'] == variant_id) & (variantid_identity_lineage_df['identity'] >= identity))].copy()
            ###########
            #
            # If some hits at this identity level, enter analysis
            #
            ###########
            if tax_lineage_by_variant_id_df.shape[0] > 0:  # If some hits at this identity, enter analysis
                ###########
                #
                # Carry out analysis if one of thise case
                # Case 1: identity >= identity_threshold
                # Case 2: identity < identity_threshold and target_tax_id.unique.count > min_number_of_taxa
                #
                ###########
                if (identity < identity_threshold and len(tax_lineage_by_variant_id_df.target_tax_id.unique().tolist()) >= min_number_of_taxa) or (identity >= identity_threshold):
                    # sort columns of tax_lineage_by_variant_id_df based on rank_hierarchy order
                    #if identity < identity_threshold:
                    #    import pdb; pdb.set_trace()
                    ltg_tax_id, ltg_rank = None, None
                    lineage_list_df_columns_sorted = [value for value in tax_lineage_by_variant_id_df if
                                                      value in rank_hierarchy]
                    tax_lineage_by_variant_id_df = tax_lineage_by_variant_id_df[lineage_list_df_columns_sorted]
                    # drop column with NaN
                    tax_lineage_by_variant_id_df = tax_lineage_by_variant_id_df.dropna(axis='columns', how='all')
                    #
                    ltg_tax_id, ltg_rank = f06_select_ltg(tax_lineage_by_variant_id_df, include_prop=include_prop)
                    # ltg_tax_id, ltg_rank = None, None
                    # if identity >= identity_threshold:
                    #     #
                    #     ltg_tax_id, ltg_rank = f06_select_ltg(tax_lineage_by_variant_id_df, include_prop=include_prop)
                    #     #
                    ###########
                    #
                    # Case 2: based on min_number_of_taxa parameter
                    #
                    ###########
                    # else:
                    #     # if tax_lineage_by_variant_id_df.shape[0] >= min_number_of_taxa:  # More/equal rows than min_number_of_taxa
                    #     if tax_lineage_by_variant_id_df.target_tax_id.unique().count() >= min_number_of_taxa:  # More/equal rows than min_number_of_taxa
                    #         #
                    #         ltg_tax_id, ltg_rank = f06_select_ltg(tax_lineage_by_variant_id_df, include_prop=include_prop)
                    if not ltg_tax_id is None:
                        lineage_dic = {}
                        lineage_dic['variant_id'] = variant_id
                        lineage_dic['identity'] = identity
                        lineage_dic['ltg_tax_id'] = ltg_tax_id
                        lineage_dic['ltg_rank'] = ltg_rank
                        # dictionnary to list
                        list_variant_id_to_ltg.append(lineage_dic)
                        break  # Do not continue lower identities
    ltg_df = pandas.DataFrame(data=list_variant_id_to_ltg)
    return ltg_df



