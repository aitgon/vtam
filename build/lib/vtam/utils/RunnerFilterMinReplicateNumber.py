import pandas


class RunnerFilterMinReplicateNumber:

    def __init__(self, variant_read_count_df):
        """Carries out a chimera analysis"""

        self.variant_read_count_df = variant_read_count_df

    def get_variant_read_count_delete_df(self, min_replicate_number):
        """
        This filter deletes variants if present in less than min_replicate_number replicates

        This filters deletes the variant if the count of the combinaison variant i and sample j
        is low then the min_replicate_number.
        The deletion condition is: count(comb (N_ij) < min_replicate_number.

        Pseudo-algorithm of this function:

        1. Compute count(comb (N_ij)
        2. Set variant/sample/replicate for deletion if count  column is low the min_replicate_number

        Updated:
        Jan 5, 2020

        :param nijk_df: Variant read count dataframe
        :type nijk_df: pandas.DataFrame

        :param min_replicate_number: Minimal number of replicates
        :type variant_read_count_input_df: int

        :return: The output of this filter is added to the 'self.variant_read_count_filter_delete_df' with 'filter_delete'=1 or 0
        :rtype: None
        """
        #
        variant_read_count_delete_df = self.variant_read_count_df.copy()
        # replicate count
        df_grouped = self.variant_read_count_df.groupby(
            by=['run_id', 'marker_id', 'variant_id', 'sample_id']).count().reset_index()
        df_grouped = df_grouped[['run_id',
                                 'marker_id',
                                 'variant_id',
                                 'sample_id',
                                 'replicate']]  # keep columns
        df_grouped = df_grouped.rename(
            columns={'replicate': 'replicate_count'})
        #
        variant_read_count_delete_df['filter_delete'] = False
        variant_read_count_delete_df = pandas.merge(
            variant_read_count_delete_df, df_grouped, on=[
                'run_id', 'marker_id', 'variant_id', 'sample_id'], how='inner')
        variant_read_count_delete_df.loc[variant_read_count_delete_df.replicate_count <
                                         min_replicate_number, 'filter_delete'] = True
        #
        variant_read_count_delete_df = variant_read_count_delete_df[[
            'run_id', 'marker_id', 'variant_id', 'sample_id', 'replicate', 'read_count', 'filter_delete']]

        return variant_read_count_delete_df
