import itertools

import pandas


class FilterRenkonenRunner(object):
    """Has attributes and methods to run the Renkonen error Filter"""

    def __init__(self, variant_read_count_df):
        """
        Initiates object for Filter Renkonen

        :param variant_read_count_df: DataFrame (run_id, marker_id, biosample_id, replicate, variant_id, read_count)
        """
        self.variant_read_count_df = variant_read_count_df

    def get_renkonen_distance(self, run_marker_biosample_df, replicate_left, replicate_right):

        """ Given run, marker, biosample and left and right replicates computes renkonen distance

        :type
        """

        replicate_pair_df = run_marker_biosample_df.loc[(run_marker_biosample_df.replicate == replicate_left) | (
                run_marker_biosample_df.replicate == replicate_right)]

        N_j_k = replicate_pair_df.groupby(by=['run_id', 'marker_id', 'biosample_id', 'replicate']).sum().reset_index()
        replicate_pair_df = replicate_pair_df.merge(N_j_k, left_on=['run_id', 'marker_id', 'biosample_id', 'replicate'],
                                                    right_on=['run_id', 'marker_id', 'biosample_id', 'replicate'])
        replicate_pair_df.rename({'read_count_x': 'N_i_j_k', 'read_count_y': 'N_j_k'}, axis=1, inplace=True)
        replicate_pair_df['N_i_j_k/N_j_k'] = replicate_pair_df['N_i_j_k'] / replicate_pair_df['N_j_k']
        replicate_pair_df = replicate_pair_df.groupby(by=['run_id', 'marker_id', 'variant_id', 'biosample_id']).min()
        renkonen_distance = 1 - replicate_pair_df['N_i_j_k/N_j_k'].sum()

        return renkonen_distance

    def get_renkonen_distance_df(self):

        renkonen_distance_df = pandas.DataFrame()

        for pos, row in self.variant_read_count_df[['run_id', 'marker_id', 'biosample_id']].drop_duplicates().iterrows():

            run_id = row.run_id
            marker_id = row.marker_id
            biosample_id = row.biosample_id

            run_marker_biosample_df = self.variant_read_count_df.loc[(self.variant_read_count_df.run_id == run_id) & (self.variant_read_count_df.marker_id == marker_id) & (
                    self.variant_read_count_df.biosample_id == biosample_id)]

            replicate_pairs = list(itertools.combinations(run_marker_biosample_df.replicate.unique(), 2))

            for replicate_left, replicate_right in replicate_pairs:

                renkonen_distance = self.get_renkonen_distance(run_marker_biosample_df,
                                                                                          replicate_left=replicate_left,
                                                                                          replicate_right=replicate_right)

                renkonen_distance_df = pandas.concat([renkonen_distance_df, pandas.DataFrame({'run_id': run_id, 'marker_id': marker_id,
                                                                  'biosample_id': biosample_id,
                                'replicate_left': replicate_left, 'replicate_right': replicate_right,
                                'renkonen_distance': renkonen_distance}, index=[0])], axis=0)

        renkonen_distance_df = renkonen_distance_df.reset_index(drop=True)

        return renkonen_distance_df

