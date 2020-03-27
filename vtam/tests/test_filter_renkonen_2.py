import itertools
import os
from unittest import TestCase

import pandas
import io

from vtam.utils.PathManager import PathManager


class TestFilterRenkonen2(TestCase):

    def setUp(self):
        pass

    def test_1(self):

        test_filter_renkonen_input_tsv_path = os.path.join(PathManager.get_test_path(), 'test_files/test_filter_renkonen_input.tsv')
        test_filter_renkonen_distance_tsv_path = os.path.join(PathManager.get_test_path(), 'test_files/test_filter_renkonen_distance.tsv')
        # with open(test_filter_renkonen_tsv_path, 'r') as fin:
        df = pandas.read_csv(test_filter_renkonen_input_tsv_path, sep="\t", header=0)

        ################################################################################################################
        #
        # Start with read count df composed of 1 run, marker and biosample, 2 replicates and any number of variants
        # Outputs a decimal between 0 and 1
        #
        ################################################################################################################

        out_df = pandas.DataFrame()

        for pos, row in df[['run_id', 'marker_id', 'biosample_id']].drop_duplicates().iterrows():

            run_id = row.run_id
            marker_id = row.marker_id
            biosample_id = row.biosample_id

            biosample_df = df.loc[(df.run_id == run_id) & (df.marker_id == marker_id) & (
                    df.biosample_id == biosample_id)]

            replicate_pairs = list(itertools.combinations(biosample_df.replicate.unique(), 2))

            for replicate_left, replicate_right in replicate_pairs:

                replicate_pair_df = biosample_df.loc[(biosample_df.replicate == replicate_left) | (
                        biosample_df.replicate == replicate_right)]

                N_j_k = replicate_pair_df.groupby(by=['run_id', 'marker_id', 'biosample_id', 'replicate']).sum().reset_index()
                replicate_pair_df = replicate_pair_df.merge(N_j_k, left_on=['run_id', 'marker_id', 'biosample_id', 'replicate'], right_on=['run_id', 'marker_id', 'biosample_id', 'replicate'])
                replicate_pair_df.rename({'read_count_x': 'N_i_j_k', 'read_count_y': 'N_j_k'}, axis=1, inplace=True)
                replicate_pair_df['N_i_j_k/N_j_k'] = replicate_pair_df['N_i_j_k'] / replicate_pair_df['N_j_k']
                replicate_pair_df = replicate_pair_df.groupby(by=['run_id', 'marker_id', 'variant_id', 'biosample_id']).min()
                renkonen_distance = 1 - replicate_pair_df['N_i_j_k/N_j_k'].sum()
                out_df = pandas.concat([out_df, pandas.DataFrame({'run_id': run_id, 'marker_id': marker_id,
                                                                  'biosample_id': biosample_id,
                                'replicate_left': replicate_left, 'replicate_right': replicate_left,
                                'renkonen_distance': renkonen_distance}, index=[0])], axis=0)
        out_df = out_df.reset_index(drop=True)

        out_bak_df = pandas.read_csv(test_filter_renkonen_distance_tsv_path, sep="\t", header=0)
        self.assertTrue(out_bak_df.to_string() == out_df.to_string())


