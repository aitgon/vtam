import os
import pandas
import unittest

from vtam.utils.RunnerFilterRenkonen import RunnerFilterRenkonen


class TestFilterRenkonen(unittest.TestCase):

    def setUp(self):
        test_filter_renkonen_input_tsv_path = os.path.join(os.path.dirname(__file__),
                                                           'test_filter_renkonen_input.tsv')
        self.variant_read_count_input_df = pandas.read_csv(test_filter_renkonen_input_tsv_path,
                                                           sep="\t", header=0)

        self.test_filter_renkonen_distance_tsv_path = os.path.join(
            os.path.dirname(__file__), 'test_filter_renkonen_distance.tsv')
        self.renkonen_distance_quantile = 0.9

        test_filter_renkonen_output_tsv_path = os.path.join(os.path.dirname(__file__),
                                                            'test_filter_renkonen_output.tsv')
        self.variant_read_count_output_df = pandas.read_csv(test_filter_renkonen_output_tsv_path,
                                                            sep="\t", header=0)

        self.filter_renkonen_runner_obj = RunnerFilterRenkonen(self.variant_read_count_input_df)

    def test_renkonen_distance(self):
        run = 1
        marker = 1
        sample = '14Ben01'
        replicate_left = '14Ben01-R1'
        replicate_right = '14Ben01-R2'

        run_marker_sample_df = self.variant_read_count_input_df.loc[
            (self.variant_read_count_input_df.run_id == run) & (
                    self.variant_read_count_input_df.marker_id == marker) & (
                    self.variant_read_count_input_df.sample_id == sample)]
        renkonen_distance = self.filter_renkonen_runner_obj.get_renkonen_distance_for_one_replicate_pair(
            run_marker_sample_df, replicate_left=replicate_left, replicate_right=replicate_right)

        self.assertTrue(renkonen_distance, 0.017426040192591086)

    def test_renkonen_distance2(self):
        run = 1
        marker = 1
        sample = '14Mon01'
        replicate_left = '14Mon01-R1'
        replicate_right = '14Mon01-R3'

        run_marker_sample_df = self.variant_read_count_input_df.loc[
            (self.variant_read_count_input_df.run_id == run) & (
                    self.variant_read_count_input_df.marker_id == marker) & (
                    self.variant_read_count_input_df.sample_id == sample)]
        renkonen_distance = self.filter_renkonen_runner_obj.get_renkonen_distance_for_one_replicate_pair(
            run_marker_sample_df,
            replicate_left=replicate_left, replicate_right=replicate_right)
        self.assertTrue(renkonen_distance, 0.375254885809898)

    def test_renkonen_distance_df(self):
        out_bak_df = pandas.read_csv(self.test_filter_renkonen_distance_tsv_path, sep="\t",
                                     header=0)
        renkonen_distance_df = self.filter_renkonen_runner_obj.get_renkonen_distance_df_for_all_sample_replicates()
        self.assertTrue(out_bak_df.to_string() == renkonen_distance_df.to_string())

    def test_get_filter_output_df(self):
        filter_out_df = self.filter_renkonen_runner_obj.get_variant_read_count_delete_df(
            self.renkonen_distance_quantile)

        filter_out_without_filter_delete_column_df = filter_out_df.loc[
            ~filter_out_df.filter_delete, ['run_id', 'marker_id', 'variant_id', 'sample_id',
                                           'replicate',
                                           'read_count']]

        self.assertTrue(self.variant_read_count_output_df.to_string(index=False)
                        == filter_out_without_filter_delete_column_df.to_string(index=False))
