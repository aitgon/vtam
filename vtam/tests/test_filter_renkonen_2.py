import itertools
import os
from unittest import TestCase

import pandas

from vtam.utils.FilterRenkonenRunner import FilterRenkonenRunner
from vtam.utils.PathManager import PathManager


class TestFilterRenkonen2(TestCase):

    @classmethod
    def setUpClass(cls):
        test_filter_renkonen_input_tsv_path = os.path.join(PathManager.get_test_path(), 'test_files/test_filter_renkonen_input.tsv')
        cls.test_filter_renkonen_distance_tsv_path = os.path.join(PathManager.get_test_path(), 'test_files/test_filter_renkonen_distance.tsv')
        cls.upper_renkonen_tail = 0.1

        cls.variant_read_count_df = pandas.read_csv(test_filter_renkonen_input_tsv_path, sep="\t", header=0)
        cls.filter_renkonen_runner_obj = FilterRenkonenRunner(cls.variant_read_count_df)
        cls.renkonen_distance_df = cls.filter_renkonen_runner_obj.get_renkonen_distance_df()

    def test_renkonen_distance(self):

        run = 1
        marker = 1
        biosample = '14Ben01'
        run_marker_biosample_df = self.variant_read_count_df.loc[
            (self.variant_read_count_df.run_id == run) & (self.variant_read_count_df.marker_id == marker) & (
                    self.variant_read_count_df.biosample_id == biosample)]
        renkonen_distance = self.filter_renkonen_runner_obj.get_renkonen_distance(run_marker_biosample_df,
                                                              replicate_left='14Ben01-R1', replicate_right='14Ben01-R2')

        self.assertTrue(renkonen_distance, 0.017426040192591086)

    def test_renkonen_distance_df(self):

        out_bak_df = pandas.read_csv(self.test_filter_renkonen_distance_tsv_path, sep="\t", header=0)
        # import pdb; pdb.set_trace()
        self.assertTrue(out_bak_df.to_string() == self.renkonen_distance_df.to_string())

    def test_1(self):

        filter_out_df = self.variant_read_count_df.copy()
        filter_out_df['filter_delete'] = False

        nb_of_replicates_df = self.variant_read_count_df[['run_id', 'marker_id', 'biosample_id', 'replicate']].drop_duplicates().groupby(
            ['run_id', 'marker_id', 'biosample_id']).count().reset_index()
        nb_of_replicates_df.rename({'replicate': 'nb_replicates'}, axis=1, inplace=True)

        self.renkonen_distance_df[
            'above_upper_renkonen_tail'] = self.renkonen_distance_df.renkonen_distance > self.upper_renkonen_tail
        for run_marker_biosample_replicate_row_i, run_marker_biosample_replicate_row in \
                self.variant_read_count_df[['run_id', 'marker_id', 'biosample_id', 'replicate']].drop_duplicates().iterrows():
            run_id = run_marker_biosample_replicate_row.run_id
            marker_id = run_marker_biosample_replicate_row.marker_id
            biosample_id = run_marker_biosample_replicate_row.biosample_id
            replicate = run_marker_biosample_replicate_row.replicate

            run_marker_biosample_replicate_renkonen_df = self.renkonen_distance_df.loc[
                (self.renkonen_distance_df.run_id == run_id) &
                (self.renkonen_distance_df.marker_id == marker_id) &
                (self.renkonen_distance_df.biosample_id == biosample_id) &
                ((self.renkonen_distance_df.replicate_left == replicate) | (self.renkonen_distance_df.replicate_right == replicate))
            ]

            nb_replicates = int(nb_of_replicates_df.loc[
                (nb_of_replicates_df.run_id == run_id) & (nb_of_replicates_df.marker_id == marker_id) & (
                            nb_of_replicates_df.biosample_id == biosample_id), 'nb_replicates'].values)

            delete_replicate = run_marker_biosample_replicate_renkonen_df.above_upper_renkonen_tail.sum() > (nb_replicates - 1) / 2
            filter_out_df.loc[(filter_out_df.run_id == run_id) & (filter_out_df.marker_id == marker_id) & (
                            filter_out_df.biosample_id == biosample_id) & (filter_out_df.replicate == replicate), 'filter_delete'] = delete_replicate

        # run_marker_biosample_row_df = self.renkonen_distance_df.merge(pandas.DataFrame(run_marker_biosample_row).T,
        #                                                               on=['run_id', 'marker_id', 'biosample_id'])
        import pdb; pdb.set_trace()
