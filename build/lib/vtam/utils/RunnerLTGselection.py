import pandas

from vtam.utils.FileParams import FileParams
from vtam.utils.constants import rank_hierarchy, identity_list


class RunnerLTGselection(object):
    """Takes a DF with columns: variant_id, %identity, 'target_tax_id', 'no rank', 'species', ...
    and the returns the LTG for each variant"""

    def __init__(self, variant_identity_lineage_df, taxonomy_df, params):

        self.variantid_identity_lineage_df = variant_identity_lineage_df
        self.taxonomy_df = taxonomy_df

        #######################################################################
        #
        # Parameters
        #
        #######################################################################

        params_dic = FileParams(params).get_params_dic()

        self.ltg_rule_threshold = params_dic['ltg_rule_threshold']
        self.include_prop = params_dic['include_prop']
        self.min_number_of_taxa = params_dic['min_number_of_taxa']

    def one_variant_to_ltg(self, variant_id):
        """Takes a DataFrame with columns: blast target identity and rank columns (lineage) for a single variant_id
        and returns a dict with values 'identity', 'ltg_tax_id', 'ltg_tax_name', 'ltg_rank']
        or an empty dictionnary if no results

        Example of the input variant_identity_lineage_df
        variant_id  identity  target_tax_id  class  family   genus  infraclass  kingdom  no rank  order  phylum  species  subclass  suborder  subphylum  superkingdom
    0           3   100.000         189839  50557  172515  189838       33339    33208   131567  30073    6656   189839      7496    185809       6960          2759
    1           3    99.429         189839  50557  172515  189838       33339    33208   131567  30073    6656   189839      7496    185809       6960          2759
    2           3    99.429         189839  50557  172515  189838       33339    33208   131567  30073    6656   189839      7496    185809       6960          2759
    3           3    99.429         189839  50557  172515  189838       33339    33208   131567  30073    6656   189839      7496    185809       6960          2759
    4           3    99.429         189839  50557  172515  189838       33339    33208   131567  30073    6656   189839      7496    185809       6960          2759

        Example of output {'identity': 80, 'ltg_tax_id': 131567.0, 'ltg_tax_name': 'cellular organisms', 'ltg_rank': 'no rank'}

        Parameters
        ----------
        blast_lineage_df int
        Variant_id to be computed the LTG

        Returns
        -------
        dict
        Return a dict with values 'identity', 'ltg_tax_id', 'ltg_tax_name', 'ltg_rank'

        """

        blast_lineage_df = self.variantid_identity_lineage_df.loc[self.variantid_identity_lineage_df['variant_id'] == variant_id]

        ltg_dic = {}

        for identity in identity_list:  # For each variant, loop each decreasing identity
            # select hits of this variant id above identity cutoff
            blast_lineage_identity_df = blast_lineage_df.loc[self.variantid_identity_lineage_df['identity'] >= identity]

            # If no hits, continue with next identity
            if blast_lineage_identity_df.shape[0] <= 0:

                continue

            ###################################################################
            #
            # Check conditions to go for the LTG
            #
            ###################################################################

            blast_target_count = len(blast_lineage_identity_df.target_tax_id.unique().tolist())

            condition_high_similarity = identity >= self.ltg_rule_threshold
            condition_low_similarity = identity < self.ltg_rule_threshold \
                                       and blast_target_count >= self.min_number_of_taxa

            if condition_low_similarity or condition_high_similarity:

                lineage_list_df_columns_sorted = [value for value in blast_lineage_identity_df if value in rank_hierarchy]
                tax_lineage_df = (blast_lineage_identity_df[lineage_list_df_columns_sorted]).copy()

                ###############################################################
                #
                # Run include_prop LTG method
                #
                ###############################################################

                ltg_tax_id, ltg_rank = self.select_ltg_include_prop(tax_lineage_df)

                if not (ltg_tax_id is None):

                    ltg_dic['identity'] = identity
                    ltg_dic['ltg_tax_id'] = ltg_tax_id
                    ltg_dic['ltg_tax_name'] = self.taxonomy_df.loc[ltg_tax_id,].name_txt
                    ltg_dic['ltg_rank'] = ltg_rank

                    return ltg_dic
        return ltg_dic

    def several_variants_to_ltg(self):
        """
        Main function that takes blast result with variant_id, target_id, identity and tax_id and returns ltg_tax_id and ltg_rank

        Example of the output variant_read_count_input_df:
        identity ltg_rank  ltg_tax_id  variant_id
    0       100  species      189839           3
    1       100  species     1077837           7
    2        99  species     1077837           9

        Args:
            variant_identity_lineage_df (pandas.DataFrame): DF with columns: variant_id, identity, target_tax_id and lineage_columns.
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

            ltg_dic = self.one_variant_to_ltg(variant_id)
            if not (ltg_dic == {}):  # if there is a ltg
                list_variant_id_to_ltg.append({**{'variant_id': variant_id}, **ltg_dic})

        ltg_df = pandas.DataFrame(data=list_variant_id_to_ltg, columns=['variant_id', 'identity',
                'ltg_tax_id', 'ltg_tax_name', 'ltg_rank'])

        ltg_df.ltg_tax_id = ltg_df.ltg_tax_id.astype('int')

        return ltg_df

    def select_ltg_include_prop(self, tax_lineage_df):
        """
        Selects LGT using the include_proc method.
        Take the LTG of the selected hits as the lowest taxonomic group contains <include_prop>
        proportion of the selected BLAST hits)

        Parameters
        ----------
        tax_lineage_df pandas DataFrame
        taxa and rank of different blast hits for a given variant id and identity
           no rank    species   genus     family     order     class  phylum
0   131567   741276.0  5533.0  1799696.0  231213.0  162481.0  5204.0
1   131567  1112827.0  6220.0   941271.0    6219.0    6218.0  6217.0
3   131567  1112827.0  6220.0   941271.0    6219.0    6218.0  6217.0
5   131567  1112827.0  6220.0   941271.0    6219.0    6218.0  6217.0
7   131567  1112827.0  6220.0   941271.0    6219.0    6218.0  6217.0
2   131567        NaN  6220.0   941271.0    6219.0    6218.0  6217.0
4   131567        NaN  6220.0   941271.0    6219.0    6218.0  6217.0
6   131567        NaN  6220.0   941271.0    6219.0    6218.0  6217.0
8   131567        NaN  6220.0   941271.0    6219.0    6218.0  6217.0

        Returns
        -------

        """

        """
        Given tax_lineage_df, selects the LTG

        Args:
            tax_lineage_df (pandas.DataFrame): DF where each column is a rank, rows are different target_ids and values are putative_ltg_ids
            include_prop (int): Percentage out of total selected qblast hits for Ltg to be present when identity>=ltg_rule_threshold

        Returns:
            int
                Taxonomical ID of Ltg
            str
                Rank of Ltg

        """

        # Remove column if all values=NaN
        tax_lineage_df.dropna(axis='columns', how='all', inplace=True)
        lineage_list_df_columns_sorted = list(
            filter(lambda x: x in tax_lineage_df.columns.tolist(),
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
        """(Pdb) putative_ltg_df
         putative_ltg_id  putative_ltg_count  putative_ltg_percentage
no rank         131567.0                   9               100.000000
phylum            6217.0                   8                88.888889
class             6218.0                   8                88.888889
order             6219.0                   8                88.888889
family          941271.0                   8                88.888889
genus             6220.0                   8                88.888889
species        1112827.0                   4                44.444444
"""
        if putative_ltg_df.putative_ltg_percentage.max() >= self.include_prop:

            ltg_tax_id = putative_ltg_df.loc[putative_ltg_df.putative_ltg_percentage >= self.include_prop, 'putative_ltg_id'].tail(1).values[0]
            ltg_rank = putative_ltg_df.loc[putative_ltg_df.putative_ltg_percentage >= self.include_prop, 'putative_ltg_id'].index[-1]

            if not (ltg_tax_id is None):

                return int(ltg_tax_id), str(ltg_rank)

        return None, None
