import pandas
from vtam.utils.FileParams import FileParams

from vtam.utils.constants import rank_hierarchy, identity_list


class RunnerLTGselection(object):
    """Takes a DF with columns: variant_id, %identity, 'target_tax_id', 'no rank', 'species', ...
    and the returns the LTG for each variant"""

    def __init__(self, variantid_identity_lineage_df, taxonomy_df, params):

        self.variantid_identity_lineage_df = variantid_identity_lineage_df
        self.taxonomy_df = taxonomy_df

        #######################################################################
        #
        # Parameters
        #
        #######################################################################

        # params_dic = constants.get_params_default_dic()
        params_dic = FileParams(params).get_params_dic()

        self.ltg_rule_threshold = params_dic['ltg_rule_threshold']
        self.include_prop = params_dic['include_prop']
        self.min_number_of_taxa = params_dic['min_number_of_taxa']

    def blast_output_to_ltg_tax_id(self):
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
            self.variantid_identity_lineage_df.variant_id.unique().tolist())
        #
        for variant_id in variant_id_list:  #  Loop sorted each variant
            for identity in identity_list:  # For each variant, loop each decreasing identity
                # select hits of this variant id above identity cutoff
                tax_lineage_by_variant_id_df = self.variantid_identity_lineage_df.loc[
                    ((self.variantid_identity_lineage_df['variant_id'] == variant_id)
                     & (self.variantid_identity_lineage_df['identity'] >= identity))].copy()
                ###########
                #
                # If some hits at this identity level, enter analysis
                #
                ###########
                if tax_lineage_by_variant_id_df.shape[0] > 0:
                    ###########
                    #
                    # Carry out analysis if one of these cases
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
                        tax_lineage_by_variant_id_df = tax_lineage_by_variant_id_df[lineage_list_df_columns_sorted]
                        # drop column with NaN
                        tax_lineage_by_variant_id_df = tax_lineage_by_variant_id_df.dropna(
                            axis='columns', how='all')
                        #
                        ltg_tax_id, ltg_rank = self.select_ltg(tax_lineage_by_variant_id_df)
                        # ltg_tax_id, ltg_rank = None, None
                        # if identity >= ltg_rule_threshold:
                        #     #
                        #     ltg_tax_id, ltg_rank = select_ltg(tax_lineage_by_variant_id_df, include_prop=include_prop)
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
                        #         ltg_tax_id, ltg_rank = select_ltg(tax_lineage_by_variant_id_df, include_prop=include_prop)
                        if not (ltg_tax_id is None):
                            lineage_dic = {}
                            lineage_dic['variant_id'] = variant_id
                            lineage_dic['identity'] = identity
                            lineage_dic['ltg_tax_id'] = ltg_tax_id
                            lineage_dic['ltg_tax_name'] = self.taxonomy_df.loc[ltg_tax_id, ].name_txt
                            lineage_dic['ltg_rank'] = ltg_rank

                            # dictionnary to list
                            list_variant_id_to_ltg.append(lineage_dic)
                            break  # Do not continue lower identities
        ltg_df = pandas.DataFrame(data=list_variant_id_to_ltg, columns=['variant_id', 'identity',
                'ltg_tax_id', 'ltg_tax_name', 'ltg_rank'])
        ltg_df.ltg_tax_id = ltg_df.ltg_tax_id.astype('int')
        return ltg_df

    def select_ltg(self, tax_lineage_df):
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

