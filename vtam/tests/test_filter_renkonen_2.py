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
        test_filter_renkonen_output_tsv_path = os.path.join(PathManager.get_test_path(), 'test_files/test_filter_renkonen_output.tsv')
        cls.upper_renkonen_tail = 0.1

        cls.variant_read_count_input_df = pandas.read_csv(test_filter_renkonen_input_tsv_path, sep="\t", header=0)
        cls.variant_read_count_output_df = pandas.read_csv(test_filter_renkonen_output_tsv_path, sep="\t", header=0)
        cls.filter_renkonen_runner_obj = FilterRenkonenRunner(cls.variant_read_count_input_df)
        cls.renkonen_distance_df = cls.filter_renkonen_runner_obj.get_renkonen_distance_df_for_all_biosample_replicates()

    def test_renkonen_distance(self):

        run = 1
        marker = 1
        biosample = '14Ben01'
        run_marker_biosample_df = self.variant_read_count_input_df.loc[
            (self.variant_read_count_input_df.run_id == run) & (self.variant_read_count_input_df.marker_id == marker) & (
                    self.variant_read_count_input_df.biosample_id == biosample)]
        renkonen_distance = self.filter_renkonen_runner_obj.get_renkonen_distance_for_one_replicate_pair(run_marker_biosample_df,
                                                                                                         replicate_left='14Ben01-R1', replicate_right='14Ben01-R2')

        self.assertTrue(renkonen_distance, 0.017426040192591086)

    def test_renkonen_distance_df(self):

        out_bak_df = pandas.read_csv(self.test_filter_renkonen_distance_tsv_path, sep="\t", header=0)
        # import pdb; pdb.set_trace()
        self.assertTrue(out_bak_df.to_string() == self.renkonen_distance_df.to_string())

    def test_get_filter_output_df(self):

        filter_out_df = self.filter_renkonen_runner_obj.get_filter_output_df(self.upper_renkonen_tail)

        filter_out_without_filter_delete_column_df = filter_out_df.loc[
            ~filter_out_df.filter_delete, ['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate',
                                           'read_count']]

        self.assertTrue(self.variant_read_count_output_df.to_string(index=False)
                        == filter_out_without_filter_delete_column_df.to_string(index=False))
