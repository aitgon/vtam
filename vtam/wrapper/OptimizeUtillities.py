import pandas


class KnownVariantAnalyzer(object):

    def __init__(self, variant_known_df, variant_read_count_df):
        """

        :param variant_known_df: This is the variant known df with columns
            - run_id
            - marker_id
            - biosample_id
            - variant_id
            - biosample_type: mock, negative, real
            - action: keep, tolerate, delete
        :param variant_biosample_df:
        """
        self.variant_known_df = variant_known_df
        self.variant_read_count_df = variant_read_count_df

    def get_variant_keep_df(self):
        variant_keep_df = self.variant_known_df.loc[
            (self.variant_known_df.action == 'keep') | (self.variant_known_df.action == 'tolerate')]
        variant_keep_df = variant_keep_df.merge(self.variant_known_df,
                                                                 on=['run_id', 'marker_id', 'biosample_id', 'variant_id'])
        variant_keep_df = variant_keep_df[
            ['run_id', 'marker_id', 'biosample_id', 'variant_id']].drop_duplicates(inplace=False)
        # Change types to int
        variant_keep_df.variant_id = variant_keep_df.variant_id.astype('uint32')
        return variant_keep_df


    def get_variant_delete_df(self):
        ##########################################################
        #
        # Get delete variants, that are not keep in mock samples
        #
        ##########################################################
        biosample_mock_list = self.variant_known_df.loc[
            self.variant_known_df.biosample_type == 'mock', 'biosample_id'].unique().tolist()
        variant_delete_mock_df = self.variant_read_count_df.loc[self.variant_read_count_df.biosample_id.isin(biosample_mock_list)]
        variant_delete_mock_df = variant_delete_mock_df.merge(self.variant_known_df, on=['run_id', 'marker_id',
                                                                                    'biosample_id', 'variant_id'],
                                                              how='left')
        variant_delete_mock_df = variant_delete_mock_df[variant_delete_mock_df.action.isnull()]
        variant_delete_mock_df = variant_delete_mock_df[['run_id', 'marker_id', 'biosample_id', 'variant_id']] \
            .drop_duplicates(inplace=False)

        ##########################################################
        #
        # Get delete variants, that are in negative samples
        #
        ##########################################################
        variant_delete_negative_df = self.variant_known_df.loc[self.variant_known_df.biosample_type == 'negative',
                                                          ['run_id', 'marker_id', 'biosample_id']]
        variant_delete_negative_df = variant_delete_negative_df.merge(self.variant_read_count_df, on=['run_id', 'marker_id',
                                                                                                'biosample_id'])
        variant_delete_negative_df = variant_delete_negative_df[['run_id', 'marker_id', 'biosample_id', 'variant_id']] \
            .drop_duplicates(inplace=False)

        ##########################################################
        #
        # Get delete variants, that are marked so in any (real) samples
        #
        ##########################################################
        variant_delete_real_df = self.variant_known_df.loc[self.variant_known_df.action == 'delete']
        variant_delete_real_df = variant_delete_real_df[~variant_delete_real_df.variant_id.isnull()]
        variant_delete_real_df = variant_delete_real_df[
            ['run_id', 'marker_id', 'biosample_id', 'variant_id']].drop_duplicates(inplace=False)

        ##########################################################
        #
        # Merge (Vertically) the three classes of delete variants
        #
        ##########################################################
        variant_delete_df = pandas.concat([variant_delete_mock_df, variant_delete_negative_df, variant_delete_real_df])
        variant_delete_df = variant_delete_df.drop_duplicates(inplace=False)
        variant_delete_df.variant_id = variant_delete_df.variant_id.astype(int)
        variant_delete_df = variant_delete_df.reset_index(drop=True)

        return variant_delete_mock_df, variant_delete_negative_df, variant_delete_real_df, variant_delete_df
