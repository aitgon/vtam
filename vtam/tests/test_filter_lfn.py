import itertools
import os
import unittest

import pandas

from vtam.utils.FilterLFNRunner import FilterLFNrunner
from vtam.utils.PathManager import PathManager


class TestFilterLFN(unittest.TestCase):

    def setUp(self):

        self.__testdir_path = os.path.join(PathManager.get_test_path())
        self.lfn_var_threshold_specific = os.path.join(
            PathManager.get_test_path(),
            self.__testdir_path,
            "test_files",
            "lfn_var_threshold_specific.tsv")

        #
        self.variant_df = pandas.DataFrame({
            'id': [1, 22],
            'sequence_': ["tata", "tgtg"],
        })
        self.variant_read_count_df = pandas.DataFrame({
            'run_id': [1] * 150,
            'marker_id': 150 * [1],
            'biosample_id': [1, 1, 1, 2, 2, 2] * 25,
            'replicate': [1, 2, 3] * 50,
            'variant_id': [*itertools.chain(*[[l] * 6 for l in range(1, 26)])],
            # [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, ..
            'read_count': [
                10, 5, 0, 249, 58, 185,
                68, 54, 100, 0, 0, 0,
                0, 0, 0, 258, 126, 500,
                0, 0, 0, 0, 1, 0,
                0, 0, 1, 0, 0, 0,
                1524, 1815, 789, 118, 98, 50,
                1, 0, 0, 0, 0, 0,
                0, 1, 0, 0, 0, 0,
                125, 214, 20, 1284, 1789, 1913,
                0, 1, 0, 0, 1, 0,
                15, 0, 1, 0, 0, 25,
                0, 0, 2, 598, 50, 875,
                2, 60, 12, 1, 0, 0,
                1, 0, 0, 0, 0, 2,
                0, 3, 0, 0, 5, 0,
                65, 98, 152, 2, 0, 1,
                52, 74, 85, 0, 0, 0,
                1, 0, 0, 5, 0, 8,
                5, 0, 1, 0, 0, 21,
                0, 0, 0, 524, 658, 125,
                0, 0, 0, 2, 0, 10,
                25, 58, 23, 10980, 8999, 13814,
                0, 5, 0, 0, 2, 0,
                1, 0, 1, 1, 0, 284,
                0, 2, 0, 0, 5, 0,
            ],
        })
        self.marker_id = 1
        #
        self.filter_lfn_runner = FilterLFNrunner(self.variant_read_count_df)

    def test_filter_lfn_runner_run(self):
        lfn_variant_threshold = 0.001
        # Occurrence is deleted if N_ijk/N_ik < lfn_variant_replicate_threshold
        # If this parameter is set (Not None), then the
        # lfn_variant_replicate_threshold instead of lfn_variant_threshold is
        # used
        lfn_variant_replicate_threshold = None
        # Occurrence is deleted if N_ijk/N_jk < lfn_ biosample
        # _replicate_threshold
        lfn_biosample_replicate_threshold = 0.001
        # Occurrence is deleted if N_ijk < lfn_ read_count_threshold
        lfn_read_count_threshold = 10

        filter_output_df = self.filter_lfn_runner.get_variant_read_count_delete_df(
            lfn_variant_threshold,
            lfn_variant_replicate_threshold,
            lfn_biosample_replicate_threshold,
            lfn_read_count_threshold)
        self.assertTrue(
            filter_output_df.filter_delete.tolist()[
                :12] == [
                0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1])

    def test_mark_delete_lfn_per_Ni(self):

        lfn_per_variant_threshold = 0.001

        self.filter_lfn_runner.mark_delete_lfn_per_Ni_or_Nik_or_Njk(
            lfn_denominator='N_i', threshold=lfn_per_variant_threshold)
        #
        self.assertTrue(self.filter_lfn_runner.variant_read_count_filter_delete_df.loc[
            (self.filter_lfn_runner.variant_read_count_filter_delete_df.variant_id == 22)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.biosample_id == 1)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.replicate == 1)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.filter_id == 2),
            'filter_delete'].values[0])
        self.assertTrue(not self.filter_lfn_runner.variant_read_count_filter_delete_df.loc[
            (self.filter_lfn_runner.variant_read_count_filter_delete_df.variant_id == 22)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.biosample_id == 1)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.replicate == 2)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.filter_id == 2),
            'filter_delete'].values[0])
        self.assertTrue(self.filter_lfn_runner.variant_read_count_filter_delete_df.loc[
            (self.filter_lfn_runner.variant_read_count_filter_delete_df.variant_id == 22)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.biosample_id == 1)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.replicate == 3)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.filter_id == 2),
            'filter_delete'].values[0])

    def test_mark_delete_lfn_per_Nik(self):

        self.filter_lfn_runner.mark_delete_lfn_per_Ni_or_Nik_or_Njk(
            lfn_denominator='N_ik', threshold=0.005)
        #
        self.assertTrue(self.filter_lfn_runner.variant_read_count_filter_delete_df.loc[
            (self.filter_lfn_runner.variant_read_count_filter_delete_df.variant_id == 12)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.biosample_id == 1)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.replicate == 3)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.filter_id == 3),
            'filter_delete'].values[0])
        self.assertTrue(not self.filter_lfn_runner.variant_read_count_filter_delete_df.loc[
            (self.filter_lfn_runner.variant_read_count_filter_delete_df.variant_id == 12)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.biosample_id == 2)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.replicate == 3)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.filter_id == 3),
            'filter_delete'].values[0])

    def test_mark_delete_lfn_absolute_read_count(self):

        self.filter_lfn_runner.mark_delete_lfn_absolute_read_count(10)
        #
        self.assertTrue(self.filter_lfn_runner.variant_read_count_filter_delete_df.loc[
            (self.filter_lfn_runner.variant_read_count_filter_delete_df.variant_id == 12)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.biosample_id == 1)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.replicate == 1)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.filter_id == 7),
            'filter_delete'].values[0])
        self.assertTrue(not self.filter_lfn_runner.variant_read_count_filter_delete_df.loc[
            (self.filter_lfn_runner.variant_read_count_filter_delete_df.variant_id == 12)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.biosample_id == 2)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.replicate == 3)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.filter_id == 7),
            'filter_delete'].values[0])
        self.assertTrue(self.filter_lfn_runner.variant_read_count_filter_delete_df.loc[
            (self.filter_lfn_runner.variant_read_count_filter_delete_df.variant_id == 1)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.biosample_id == 1)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.replicate == 2)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.filter_id == 7),
            'filter_delete'].values[0])

    def test_mark_delete_lfn_per_Njk(self):

        self.filter_lfn_runner.mark_delete_lfn_per_Ni_or_Nik_or_Njk(
            lfn_denominator='N_jk', threshold=0.001)
        #
        self.assertTrue(
            not self.filter_lfn_runner.variant_read_count_filter_delete_df.loc[
                (self.filter_lfn_runner.variant_read_count_filter_delete_df.variant_id == 9) & (
                    self.filter_lfn_runner.variant_read_count_filter_delete_df.biosample_id == 2) & (
                    self.filter_lfn_runner.variant_read_count_filter_delete_df.replicate == 3) & (
                    self.filter_lfn_runner.variant_read_count_filter_delete_df.filter_id == 6),
                'filter_delete'].values[0])
        self.assertTrue(self.filter_lfn_runner.variant_read_count_filter_delete_df.loc[
            (self.filter_lfn_runner.variant_read_count_filter_delete_df.variant_id == 12)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.biosample_id == 1)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.replicate == 1)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.filter_id == 6),
            'filter_delete'].values[0])
        self.assertTrue(not self.filter_lfn_runner.variant_read_count_filter_delete_df.loc[
            (self.filter_lfn_runner.variant_read_count_filter_delete_df.variant_id == 12)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.biosample_id == 1)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.replicate == 3)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.filter_id == 6),
            'filter_delete'].values[0])

        self.assertTrue(self.filter_lfn_runner.variant_read_count_filter_delete_df.loc[
            (self.filter_lfn_runner.variant_read_count_filter_delete_df.variant_id == 24)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.biosample_id == 1)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.replicate == 3)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.filter_id == 6),
            'filter_delete'].values[0])

        self.assertTrue(self.filter_lfn_runner.variant_read_count_filter_delete_df.loc[
            (self.filter_lfn_runner.variant_read_count_filter_delete_df.variant_id == 24)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.biosample_id == 1)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.replicate == 1)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.filter_id == 6),
            'filter_delete'].values[0])

    def test_mark_delete_lfn_do_not_pass_all_filters(self):
        self.filter_lfn_runner.mark_delete_lfn_per_Ni_or_Nik_or_Njk(
            lfn_denominator='N_i', threshold=0.001)

        self.filter_lfn_runner.mark_delete_lfn_per_Ni_or_Nik_or_Njk(
            lfn_denominator='N_ik', threshold=0.005)

        self.filter_lfn_runner.mark_delete_lfn_do_not_pass_all_filters()

        self.assertTrue(not self.filter_lfn_runner.variant_read_count_filter_delete_df.loc[
            (self.filter_lfn_runner.variant_read_count_filter_delete_df.run_id == 1)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.variant_id == 1)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.biosample_id == 1)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.replicate == 1)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.filter_id == 8),
            'filter_delete'].values[0])
        self.assertTrue(self.filter_lfn_runner.variant_read_count_filter_delete_df.loc[
            (self.filter_lfn_runner.variant_read_count_filter_delete_df.run_id == 1) & (
                self.filter_lfn_runner.variant_read_count_filter_delete_df.variant_id == 1)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.biosample_id == 1)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.replicate == 3)
            & (self.filter_lfn_runner.variant_read_count_filter_delete_df.filter_id == 8),
            'filter_delete'].values[0])
