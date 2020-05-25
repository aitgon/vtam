import inspect
import os
import pathlib

import pandas
from Bio.Blast.Applications import NcbiblastnCommandline

from vtam.utils.Logger import Logger
from vtam.utils.PathManager import PathManager
from vtam.utils.constants import rank_hierarchy, identity_list


class TaxAssignRunner(object):
    """Will assign variants to a taxon"""

    def __init__(
            self,
            sequence_list,
            taxonomy_df,
            blast_db_dir,
            blast_db_name,
            ltg_rule_threshold,
            include_prop,
            min_number_of_taxa,
            num_threads,
            qcov_hsp_perc):

        # self.variant_df = variant_df
        # stores tax_id and old_tax_id
        self.old_tax_id_df = taxonomy_df[['old_tax_id']].drop_duplicates()
        self.taxonomy_df = taxonomy_df[[
            'parent_tax_id', 'rank', 'name_txt']].drop_duplicates()
        self.blast_db_dir = blast_db_dir
        self.this_temp_dir = os.path.join(
            PathManager.instance().get_tempdir(),
            os.path.basename(__file__))
        pathlib.Path(self.this_temp_dir).mkdir(exist_ok=True)
        self.ltg_rule_threshold = ltg_rule_threshold
        self.include_prop = include_prop
        self.min_number_of_taxa = min_number_of_taxa
        self.num_threads = num_threads
        self.qcov_hsp_perc = qcov_hsp_perc

        #######################################################################
        #
        # 2 Create FASTA file with Variants
        #
        #######################################################################

        Logger.instance().debug(
            "file: {}; line: {}; Create SortedReadFile from Variants".format(
                __file__, inspect.currentframe().f_lineno))
        variant_fasta = os.path.join(self.this_temp_dir, 'variant.fasta')
        # variant_df_utils = VariantDF(variant_df)
        # variant_df_utils.to_fasta(variant_fasta)
        with open(variant_fasta, 'w') as fout:
            for seq in sequence_list:
                fout.write(">{}\n{}\n".format(seq, seq))

        #######################################################################
        #
        # 3 Run local blast
        #
        #######################################################################

        Logger.instance().debug(
            "file: {}; line: {}; Running local blast with FASTA input {}".format(
                __file__, inspect.currentframe().f_lineno, variant_fasta))

        # Run and read local blast result
        blast_output_tsv = os.path.join(self.this_temp_dir, 'blast_output.tsv')
        # blast_output_tsv = "/home/gonzalez/tmp/blast/blast_output.tsv" # uncomment for testing
        # get blast db dir and filename prefix from NHR file
        os.environ['BLASTDB'] = self.blast_db_dir

        blastn_cline = NcbiblastnCommandline(
            query=variant_fasta,
            db=blast_db_name,
            evalue=1e-5,
            outfmt='"6 qseqid sacc pident evalue qcovhsp staxids"',
            dust='yes',
            qcov_hsp_perc=self.qcov_hsp_perc,
            num_threads=self.num_threads,
            out=blast_output_tsv)
        Logger.instance().debug(
            "file: {}; line: {}; {}".format(
                __file__,
                inspect.currentframe().f_lineno,
                str(blastn_cline)))
        #
        # Run blast
        stdout, stderr = blastn_cline()

        #######################################################################
        #
        # Process blast results
        #
        #######################################################################

        Logger.instance().debug(
            "file: {}; line: {}; Reading Blast output from: {}".format(
                __file__, inspect.currentframe().f_lineno, blast_output_tsv))
        blast_output_df = pandas.read_csv(
            blast_output_tsv,
            sep='\t',
            header=None,
            names=[
                'variant_id',
                'target_id',
                'identity',
                'evalue',
                'coverage',
                'target_tax_id'])
        # Remove null target tax ids
        blast_output_df = blast_output_df.loc[~blast_output_df.target_tax_id.isnull(
        )]
        # expand multiple target_tax_ids
        blast_output_df.target_tax_id = blast_output_df.target_tax_id.astype(
            'str')  # first convert as string
        blast_output_df = (
            pandas.concat(
                [
                    blast_output_df,
                    blast_output_df.target_tax_id.str.split(
                        pat=';',
                        n=1,
                        expand=True)],
                axis=1))
        self.ltg_df = None  # Init it
        if blast_output_df.shape[0] > 0:  # If blast output exists
            # Select first tax_id
            blast_output_df = blast_output_df[[
                'variant_id', 'target_id', 'identity', 'evalue', 'coverage', 0]]
            # rename first tax_id
            blast_output_df = blast_output_df.rename(
                columns={0: 'target_tax_id'})
            # Convert columns back to int
            blast_output_df.target_tax_id = blast_output_df.target_tax_id.astype(
                'float')
            blast_output_df.target_tax_id = blast_output_df.target_tax_id.astype(
                'int')
            # Blast output extract
            """   variant_id  target_id  identity        evalue  coverage  target_tax_id
    0           2  MF7836761    99.429  1.620000e-86       100        1469487
    1           2  MF7836761    99.429  1.620000e-86       100         189839
    2           2  KY2618191    98.857  7.520000e-85       100         189839
    3           2  MF7834791    98.857  7.520000e-85       100         189839
    4           2  KU9559321    98.857  7.520000e-85       100         189839
    """

            ##########################################################
            #
            # Read target_tax_id
            # Compute lineages for each unique target_tax_id
            # Create a DF with these columns: tax_id and its lineage in wide format
            # Merge to the blast result
            #
            ##########################################################

            Logger.instance().debug(
                "file: {}; line: {}; Open taxonomy.tsv DB".format(
                    __file__, inspect.currentframe().f_lineno))
            blast_output_df.target_tax_id = pandas.to_numeric(
                blast_output_df.target_tax_id)
            # getting the taxonomy_db to variant_read_count_input_df
            # taxonomy_tsv_path = taxonomy_tsv
            #
            Logger.instance().debug(
                "file: {}; line: {}; Annotate each target_tax_id with its lineage as columns in wide format".format(
                    __file__, inspect.currentframe().f_lineno))
            lineage_list = []
            for target_tax_id_i, target_tax_id in enumerate(
                    blast_output_df.target_tax_id.unique().tolist()):
                if target_tax_id_i % 100 == 0:
                    Logger.instance().debug(
                        "Get lineage of {}-th tax id {} (Total {} tax ids)".format(
                            target_tax_id_i, target_tax_id, len(
                                blast_output_df.target_tax_id.unique().tolist())))
                lineage_list.append(
                    self.tax_id_to_taxonomy_lineage(
                        tax_id=target_tax_id))
            tax_id_to_lineage_df = pandas.DataFrame(lineage_list)
            #
            Logger.instance().debug(
                "file: {}; line: {}; Merge blast result including tax_id with their lineages".format(
                    __file__, inspect.currentframe().f_lineno))
            # Merge local blast output with tax_id_to_lineage_df
            variantid_identity_lineage_df = blast_output_df.merge(
                tax_id_to_lineage_df, left_on='target_tax_id', right_on='tax_id')
            variantid_identity_lineage_df.drop('tax_id', axis=1, inplace=True)

            ##########################################################
            #
            #  6 test_f05_select_ltg_identity
            #
            ##########################################################

            Logger.instance().debug(
                "file: {}; line: {}; Main loop over variant and identity to"
                "compute the whole set of ltg_tax_id and ltg_rank for each variant_id"
                "to a dataframe".format(
                    __file__, inspect.currentframe().f_lineno))
            #
            # f07_blast_result_to_ltg_tax_id(tax_lineage_df,ltg_rule_threshold=97, include_prop=90, min_number_of_taxa=3):
            # this function return a data frame containing the Ltg rank and Ltg Tax_id for each variant
            #
            self.ltg_df = self.f07_blast_result_to_ltg_tax_id(
                variantid_identity_lineage_df)

    def tax_id_to_taxonomy_lineage(self, tax_id):
        """
        Takes a tax_id and creates a dictionary with the taxonomy lineage

        :param tax_id: Identifier of taxon
        :type tax_id: int
        :param return_tax_name: (pandas.DataFrame: DataFrame with taxonomy information
        :type return_tax_name: str

        :return: Dic with ranks/keys and tax_id/values
        :rtype: list
        Args:
            tax_id (Integer): Identifier of taxon
            taxonomy_df (pandas.DataFrame: DataFrame with taxonomy information

        Returns:
            Dictionnary: with taxonomy information for given tax_id, {'tax_id': 183142, 'species': 183142, 'genus': 10194, 'family': 10193, 'order': 84394,
                             'superorder': 1709201, 'class': 10191, 'phylum': 10190, 'no rank': 131567, 'kingdom': 33208,
                             'superkingdom': 2759}

        """

        lineage_dic = {'tax_id': tax_id}
        # while tax_id != 1:
        #     # try to use taxonomy_df.tax_id
        #     # import pdb; pdb.set_trace()
        #     # tax_id_row = self.taxonomy_df.loc[self.taxonomy_df.tax_id == tax_id, ]
        #     if tax_id in self.taxonomy_df.index:
        #
        #     try:
        #         tax_id_row = self.taxonomy_df.loc[tax_id, ]
        #     except KeyError:
        #     import pdb; pdb.set_trace()
        #     # row empty, try to use old_tax_id
        #     if tax_id_row.shape[0] == 0:
        #         tax_id_row = self.taxonomy_df.loc[self.taxonomy_df.old_tax_id == tax_id, ]
        #     rank = tax_id_row['rank'].values[0]
        #     parent_tax_id = tax_id_row['parent_tax_id'].values[0]
        #     lineage_dic[rank] = tax_id
        #     # if return_tax_name: # return tax_name instead of tax_id
        #     #     tax_name = tax_id_row['name_txt'].values[0]
        #     #     lineage_dic[rank] = tax_name
        #     tax_id = parent_tax_id
        while tax_id != 1:
            if tax_id in self.taxonomy_df.index:
                tax_id_row = self.taxonomy_df.loc[tax_id, ]
            elif tax_id in self.old_tax_id_df.old_tax_id.tolist():  # Try old tax id
                tax_id2 = self.old_tax_id_df.loc[self.old_tax_id_df.old_tax_id ==
                                                 tax_id, 'old_tax_id'].index[0]
                tax_id_row = self.taxonomy_df.loc[tax_id2, ]
            else:
                return lineage_dic
            rank = tax_id_row['rank']
            parent_tax_id = tax_id_row['parent_tax_id']
            lineage_dic[rank] = tax_id
            tax_id = parent_tax_id
        return lineage_dic

    def f06_select_ltg(self, tax_lineage_df):
        """
        Given tax_lineage_df, selects the LTG

        Args:
            tax_lineage_df (pandas.DataFrame): DF where each column is a rank, rows are different target_ids and values are putative_ltg_ids
            include_prop (int): Percentage out of total selected qblast hits for Ltg to be present when identity>=ltg_rule_threshold

        Returns:
            ltg_tax_id (int): Taxonomical ID of Ltg
            ltg_rank (str): Rank of Ltg

        """

        # Remove column if all values=NaN
        tax_lineage_df.dropna(axis='columns', how='all', inplace=True)
        lineage_list_df_columns_sorted = list(
            filter(
                lambda x: x in tax_lineage_df.columns.tolist(),
                rank_hierarchy))
        tax_lineage_df = tax_lineage_df[lineage_list_df_columns_sorted]
        putative_ltg_df = pandas.DataFrame(
            {
                'putative_ltg_id': tax_lineage_df.apply(
                    lambda x: x.value_counts().index[0],
                    axis=0),
                'putative_ltg_count': tax_lineage_df.apply(
                    lambda x: x.value_counts().iloc[0],
                    axis=0)})
        putative_ltg_df['putative_ltg_percentage'] = putative_ltg_df.putative_ltg_count / \
            tax_lineage_df.shape[0] * 100
        ltg_tax_id = putative_ltg_df.loc[putative_ltg_df.putative_ltg_percentage >=
                                         self.include_prop, 'putative_ltg_id'].tail(1).values[0]
        ltg_rank = putative_ltg_df.loc[putative_ltg_df.putative_ltg_percentage >=
                                       self.include_prop, 'putative_ltg_id'].index[-1]
        return ltg_tax_id, ltg_rank

    def f07_blast_result_to_ltg_tax_id(self, variantid_identity_lineage_df):
        """
        Main function that takes blast result with variant_id, target_id, identity and tax_id and returns ltg_tax_id and ltg_rank

        Example of the input variantid_identity_lineage_df
        variant_id  identity  target_tax_id  class  family   genus  infraclass  kingdom  no rank  order  phylum  species  subclass  suborder  subphylum  superkingdom
    0           3   100.000         189839  50557  172515  189838       33339    33208   131567  30073    6656   189839      7496    185809       6960          2759
    1           3    99.429         189839  50557  172515  189838       33339    33208   131567  30073    6656   189839      7496    185809       6960          2759
    2           3    99.429         189839  50557  172515  189838       33339    33208   131567  30073    6656   189839      7496    185809       6960          2759
    3           3    99.429         189839  50557  172515  189838       33339    33208   131567  30073    6656   189839      7496    185809       6960          2759
    4           3    99.429         189839  50557  172515  189838       33339    33208   131567  30073    6656   189839      7496    185809       6960          2759
        Example of the output variant_read_count_input_df:
        identity ltg_rank  ltg_tax_id  variant_id
    0       100  species      189839           3
    1       100  species     1077837           7
    2        99  species     1077837           9

        Args:
            variantid_identity_lineage_df (pandas.DataFrame): DF with columns: variant_id, identity, target_tax_id and lineage_columns.
            ltg_rule_threshold (int): Identity value where we change of using include_prop method to min_number_of_taxa, default 97
            include_prop (int): Percentage out of total selected qblast hits for Ltg to be present when identity>=ltg_rule_threshold
            min_number_of_taxa (int): Minimal number of taxa, where LTF must be present when identity<ltg_rule_threshold

        Returns:
            ltg_df (pandas.DataFrame): DF with variant_id, ltg_tax_id and ltg_tax_rank

        """
        #
        list_variant_id_to_ltg = []
        variant_id_list = sorted(
            variantid_identity_lineage_df.variant_id.unique().tolist())
        #
        for variant_id in variant_id_list:  #  Loop sorted each variant
            for identity in identity_list:  # For each variant, loop each decreasing identity
                tax_lineage_by_variant_id_df = variantid_identity_lineage_df.loc[
                    ((variantid_identity_lineage_df['variant_id'] == variant_id)
                     & (variantid_identity_lineage_df['identity'] >= identity))].copy()
                ###########
                #
                # If some hits at this identity level, enter analysis
                #
                ###########
                # If some hits at this identity, enter analysis
                if tax_lineage_by_variant_id_df.shape[0] > 0:
                    ###########
                    #
                    # Carry out analysis if one of thise case
                    # Case 1: identity >= ltg_rule_threshold
                    # Case 2: identity < ltg_rule_threshold and target_tax_id.unique.count > min_number_of_taxa
                    #
                    ###########
                    if (identity < self.ltg_rule_threshold and len(tax_lineage_by_variant_id_df.target_tax_id.unique(
                    ).tolist()) >= self.min_number_of_taxa) or (identity >= self.ltg_rule_threshold):
                        # sort columns of tax_lineage_by_variant_id_df based on
                        # rank_hierarchy order
                        ltg_tax_id, ltg_rank = None, None
                        lineage_list_df_columns_sorted = [
                            value for value in tax_lineage_by_variant_id_df if value in rank_hierarchy]
                        tax_lineage_by_variant_id_df = tax_lineage_by_variant_id_df[
                            lineage_list_df_columns_sorted]
                        # drop column with NaN
                        tax_lineage_by_variant_id_df = tax_lineage_by_variant_id_df.dropna(
                            axis='columns', how='all')
                        #
                        ltg_tax_id, ltg_rank = self.f06_select_ltg(
                            tax_lineage_by_variant_id_df)
                        # ltg_tax_id, ltg_rank = None, None
                        # if identity >= ltg_rule_threshold:
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
                        if ltg_tax_id is not None:
                            lineage_dic = {}
                            lineage_dic['variant_id'] = variant_id
                            lineage_dic['identity'] = identity
                            lineage_dic['ltg_tax_id'] = ltg_tax_id
                            lineage_dic['ltg_tax_name'] = self.taxonomy_df.loc[ltg_tax_id, ].name_txt
                            lineage_dic['ltg_rank'] = ltg_rank

                            # dictionnary to list
                            list_variant_id_to_ltg.append(lineage_dic)
                            break  # Do not continue lower identities
        ltg_df = pandas.DataFrame(
            data=list_variant_id_to_ltg,
            columns=[
                'variant_id',
                'identity',
                'ltg_tax_id',
                'ltg_tax_name',
                'ltg_rank'])
        return ltg_df


def f07_blast_result_to_ltg_tax_id(
        variantid_identity_lineage_df,
        ltg_rule_threshold,
        include_prop,
        min_number_of_taxa):
    """
    Main function that takes blast result with variant_id, target_id, identity and tax_id and returns ltg_tax_id and ltg_rank

    Example of the input variantid_identity_lineage_df
    variant_id  identity  target_tax_id  class  family   genus  infraclass  kingdom  no rank  order  phylum  species  subclass  suborder  subphylum  superkingdom
0           3   100.000         189839  50557  172515  189838       33339    33208   131567  30073    6656   189839      7496    185809       6960          2759
1           3    99.429         189839  50557  172515  189838       33339    33208   131567  30073    6656   189839      7496    185809       6960          2759
2           3    99.429         189839  50557  172515  189838       33339    33208   131567  30073    6656   189839      7496    185809       6960          2759
3           3    99.429         189839  50557  172515  189838       33339    33208   131567  30073    6656   189839      7496    185809       6960          2759
4           3    99.429         189839  50557  172515  189838       33339    33208   131567  30073    6656   189839      7496    185809       6960          2759
    Example of the output variant_read_count_input_df:
    identity ltg_rank  ltg_tax_id  variant_id
0       100  species      189839           3
1       100  species     1077837           7
2        99  species     1077837           9

    Args:
        variantid_identity_lineage_df (pandas.DataFrame): DF with columns: variant_id, identity, target_tax_id and lineage_columns.
        ltg_rule_threshold (int): Identity value where we change of using include_prop method to min_number_of_taxa, default 97
        include_prop (int): Percentage out of total selected qblast hits for Ltg to be present when identity>=ltg_rule_threshold
        min_number_of_taxa (int): Minimal number of taxa, where LTF must be present when identity<ltg_rule_threshold

    Returns:
        ltg_df (pandas.DataFrame): DF with variant_id, ltg_tax_id and ltg_tax_rank

    """
    #
    list_variant_id_to_ltg = []
    variant_id_list = sorted(
        variantid_identity_lineage_df.variant_id.unique().tolist())
    #
    for variant_id in variant_id_list:  #  Loop sorted each variant
        for identity in identity_list:  # For each variant, loop each decreasing identity
            tax_lineage_by_variant_id_df = variantid_identity_lineage_df.loc[(
                (variantid_identity_lineage_df['variant_id'] == variant_id) & (variantid_identity_lineage_df['identity'] >= identity))].copy()
            ###########
            #
            # If some hits at this identity level, enter analysis
            #
            ###########
            # If some hits at this identity, enter analysis
            if tax_lineage_by_variant_id_df.shape[0] > 0:
                ###########
                #
                # Carry out analysis if one of thise case
                # Case 1: identity >= ltg_rule_threshold
                # Case 2: identity < ltg_rule_threshold and target_tax_id.unique.count > min_number_of_taxa
                #
                ###########
                if (identity < ltg_rule_threshold and len(tax_lineage_by_variant_id_df.target_tax_id.unique(
                ).tolist()) >= min_number_of_taxa) or (identity >= ltg_rule_threshold):
                    # sort columns of tax_lineage_by_variant_id_df based on
                    # rank_hierarchy order
                    ltg_tax_id, ltg_rank = None, None
                    lineage_list_df_columns_sorted = [
                        value for value in tax_lineage_by_variant_id_df if value in rank_hierarchy]
                    tax_lineage_by_variant_id_df = tax_lineage_by_variant_id_df[
                        lineage_list_df_columns_sorted]
                    # drop column with NaN
                    tax_lineage_by_variant_id_df = tax_lineage_by_variant_id_df.dropna(
                        axis='columns', how='all')
                    #
                    ltg_tax_id, ltg_rank = f06_select_ltg(
                        tax_lineage_by_variant_id_df, include_prop=include_prop)
                    # ltg_tax_id, ltg_rank = None, None
                    # if identity >= ltg_rule_threshold:
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
                    if ltg_tax_id is not None:
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


def f06_select_ltg(tax_lineage_df, include_prop):
    """
    Given tax_lineage_df, selects the Ltg

    Args:
        tax_lineage_df (pandas.DataFrame): DF where each column is a rank, rows are different target_ids and values are putative_ltg_ids
        include_prop (int): Percentage out of total selected qblast hits for Ltg to be present when identity>=ltg_rule_threshold

    Returns:
        ltg_tax_id (int): Taxonomical ID of Ltg
        ltg_rank (str): Rank of Ltg

    """

    # Remove column if all values=NaN
    tax_lineage_df.dropna(axis='columns', how='all', inplace=True)
    lineage_list_df_columns_sorted = list(
        filter(
            lambda x: x in tax_lineage_df.columns.tolist(),
            rank_hierarchy))
    tax_lineage_df = tax_lineage_df[lineage_list_df_columns_sorted]
    putative_ltg_df = pandas.DataFrame(
        {
            'putative_ltg_id': tax_lineage_df.apply(
                lambda x: x.value_counts().index[0],
                axis=0),
            'putative_ltg_count': tax_lineage_df.apply(
                lambda x: x.value_counts().iloc[0],
                axis=0)})
    putative_ltg_df['putative_ltg_percentage'] = putative_ltg_df.putative_ltg_count / \
        tax_lineage_df.shape[0] * 100
    ltg_tax_id = putative_ltg_df.loc[putative_ltg_df.putative_ltg_percentage >=
                                     include_prop, 'putative_ltg_id'].tail(1).values[0]
    ltg_rank = putative_ltg_df.loc[putative_ltg_df.putative_ltg_percentage >=
                                   include_prop, 'putative_ltg_id'].index[-1]
    return ltg_tax_id, ltg_rank
