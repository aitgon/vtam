import itertools
import pandas


class RunnerFilterRenkonen(object):
    """Has attributes and methods to run_name the Renkonen error Filter"""

    def __init__(self, variant_read_count_df):
        """
        Initiates object for Filter Renkonen

        :param variant_read_count_df: DataFrame (run_id, marker_id, sample_id, replicate, variant_id, read_count)
        """
        self.variant_read_count_df = variant_read_count_df

    def get_variant_read_count_delete_df(self, renkonen_distance_quantile):

        filter_out_df = self.variant_read_count_df.copy()
        filter_out_df['filter_delete'] = False

        nb_of_replicates_df = self.variant_read_count_df[['run_id', 'marker_id', 'sample_id', 'replicate']]\
            .drop_duplicates().groupby(
            ['run_id', 'marker_id', 'sample_id']).count().reset_index()
        nb_of_replicates_df.rename(
            {'replicate': 'nb_replicates'}, axis=1, inplace=True)

        renkonen_distance_df = self.get_renkonen_distance_df_for_all_sample_replicates()

        renkonen_distance_cutoff = renkonen_distance_df.renkonen_distance.quantile(
            renkonen_distance_quantile, interpolation='lower')

        renkonen_distance_df[
            'above_renkonen_distance_quantile'] = renkonen_distance_df.renkonen_distance > renkonen_distance_cutoff
        for run_marker_sample_replicate_row_i, run_marker_sample_replicate_row in self.variant_read_count_df[[
                'run_id', 'marker_id', 'sample_id', 'replicate']].drop_duplicates().iterrows():
            run_id = run_marker_sample_replicate_row.run_id
            marker_id = run_marker_sample_replicate_row.marker_id
            sample_id = run_marker_sample_replicate_row.sample_id
            replicate = run_marker_sample_replicate_row.replicate

            run_marker_sample_replicate_renkonen_df = renkonen_distance_df.loc[
                (renkonen_distance_df.run_id == run_id) & (
                    renkonen_distance_df.marker_id == marker_id) & (
                    renkonen_distance_df.sample_id == sample_id) & (
                    (renkonen_distance_df.replicate_left == replicate) | (
                        renkonen_distance_df.replicate_right == replicate))]

            nb_replicates = int(nb_of_replicates_df.loc[
                (nb_of_replicates_df.run_id == run_id) & (nb_of_replicates_df.marker_id == marker_id) & (
                    nb_of_replicates_df.sample_id == sample_id), 'nb_replicates'].values)

            delete_replicate = run_marker_sample_replicate_renkonen_df.above_renkonen_distance_quantile.sum(
            ) > (nb_replicates - 1) / 2
            filter_out_df.loc[(filter_out_df.run_id == run_id) & (filter_out_df.marker_id == marker_id) &
                              (filter_out_df.sample_id == sample_id) &
                              (filter_out_df.replicate == replicate), 'filter_delete'] = delete_replicate
        return filter_out_df

    def get_renkonen_distance_for_one_replicate_pair(
            self, run_marker_sample_df, replicate_left, replicate_right):
        """ Given run_name, marker_name, sample and left and right replicates computes renkonen distance

        :type
        """

        replicate_pair_vertical_df = run_marker_sample_df.loc[(run_marker_sample_df.replicate == replicate_left) | (
            run_marker_sample_df.replicate == replicate_right)]

        N_j_k_df = replicate_pair_vertical_df[['run_id',
                                               'marker_id',
                                               'sample_id',
                                               'replicate',
                                               'read_count']] .groupby(by=['run_id',
                                                                           'marker_id',
                                                                           'sample_id',
                                                                           'replicate']).sum().reset_index()
        replicate_pair_vertical_df = replicate_pair_vertical_df.merge(
            N_j_k_df, left_on=[
                'run_id', 'marker_id', 'sample_id', 'replicate'], right_on=[
                'run_id', 'marker_id', 'sample_id', 'replicate'])
        replicate_pair_vertical_df.rename(
            {'read_count_x': 'N_i_j_k', 'read_count_y': 'N_j_k'}, axis=1, inplace=True)
        replicate_pair_vertical_df['N_i_j_k/N_j_k'] = replicate_pair_vertical_df['N_i_j_k'] / \
            replicate_pair_vertical_df['N_j_k']
        replicate_left_df = replicate_pair_vertical_df.loc[
            replicate_pair_vertical_df.replicate == replicate_left]
        replicate_right_df = replicate_pair_vertical_df.loc[
            replicate_pair_vertical_df.replicate == replicate_right]
        replicate_pair_horizontal_df = replicate_left_df.merge(
            replicate_right_df, on=[
                'run_id', 'marker_id', 'sample_id', 'variant_id'], how='outer')
        replicate_pair_horizontal_df.fillna(0, inplace=True)
        replicate_pair_horizontal_df['min'] = replicate_pair_horizontal_df.apply(
            lambda x: min(x['N_i_j_k/N_j_k_x'], x['N_i_j_k/N_j_k_y']), axis=1)
        renkonen_distance = 1 - replicate_pair_horizontal_df['min'].sum()

        return renkonen_distance

    def get_renkonen_distance_df_for_all_sample_replicates(self):

        renkonen_distance_df = pandas.DataFrame()

        for pos, row in self.variant_read_count_df[[
                'run_id', 'marker_id', 'sample_id']].drop_duplicates().iterrows():

            run_id = row.run_id
            marker_id = row.marker_id
            sample_id = row.sample_id

            run_marker_sample_df = self.variant_read_count_df.loc[(self.variant_read_count_df.run_id == run_id) & (
                self.variant_read_count_df.marker_id == marker_id) & (self.variant_read_count_df.sample_id == sample_id)]

            replicate_pairs = list(
                itertools.combinations(
                    run_marker_sample_df.replicate.unique(), 2))

            for replicate_left, replicate_right in replicate_pairs:

                renkonen_distance = self.get_renkonen_distance_for_one_replicate_pair(
                    run_marker_sample_df, replicate_left=replicate_left, replicate_right=replicate_right)

                renkonen_distance_df = pandas.concat(
                    [
                        renkonen_distance_df,
                        pandas.DataFrame(
                            {
                                'run_id': run_id,
                                'marker_id': marker_id,
                                'sample_id': sample_id,
                                'replicate_left': replicate_left,
                                'replicate_right': replicate_right,
                                'renkonen_distance': renkonen_distance},
                            index=[0])],
                    axis=0)

        renkonen_distance_df = renkonen_distance_df.reset_index(drop=True)

        return renkonen_distance_df
