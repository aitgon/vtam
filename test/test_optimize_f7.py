import os
import sqlite3
from unittest import TestCase

import pandas

from wopmetabarcoding.utils.PathFinder import PathFinder
from wopmetabarcoding.utils.utilities import create_step_tmp_dir
from wopmetabarcoding.wrapper.FilterLFNutilities import FilterLFNRunner
from wopmetabarcoding.wrapper.FilterPCRError import f10_pcr_error_run_vsearch


class TestOptimizeF7(TestCase):

    def setUp(self):
        self.__testdir_path = os.path.join(PathFinder.get_module_test_path())
        self.read_count_keep_delete_path = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "test_files",
                                                  "optimize_f7", "read_coutn_keep_delete.tsv")
        self.this_step_tmp_dir = create_step_tmp_dir(__file__)
        #

    def test_01(self):


        ####
        #
        #  read_count of keep and delete
        #
        ####
        read_count_keep_delete_df = pandas.read_csv(self.read_count_keep_delete_path, sep='\t', header=0)
        #
        variant_read_count_df = read_count_keep_delete_df[['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id', 'N_ijk']]
        variant_read_count_df = variant_read_count_df.rename(columns={"N_ijk": "read_count"})
        #
        variant_keep_list = (
        read_count_keep_delete_df.loc[read_count_keep_delete_df.action == 'keep']).variant_id.unique().tolist()
        variant_delete_list = (
        read_count_keep_delete_df.loc[read_count_keep_delete_df.action == 'delete']).variant_id.unique().tolist()

        ##############
        #
        # Main loop through parameter values
        #
        ##############
        lfn_read_count_threshold_min = 10
        lfn_read_count_threshold_max = lfn_read_count_threshold_min * 1000
        lfn_read_count_threshold = lfn_read_count_threshold_min
        #
        lfn_per_sum_variant_threshold_min = 0.001 # default value
        lfn_per_sum_variant_threshold_max = lfn_per_sum_variant_threshold_min * 1000
        lfn_per_sum_variant_threshold = lfn_per_sum_variant_threshold_min
        #
        out_lfn_per_sum_variant_list = []
        while lfn_per_sum_variant_threshold <= lfn_per_sum_variant_threshold_max:
            while lfn_read_count_threshold <= lfn_read_count_threshold_max:
                #
                lfn_filter_runner = FilterLFNRunner(variant_read_count_df)
                lfn_filter_runner.f7_lfn_delete_absolute_read_count(lfn_read_count_threshold)
                variant_read_cound_remained_df = lfn_filter_runner.delete_variant_df
                variant_read_cound_remained_df = variant_read_cound_remained_df.loc[
                    variant_read_cound_remained_df.filter_delete == 0]
                #
                variant_remained_list = variant_read_cound_remained_df.variant_id.unique().tolist()
                #  Count how many from 'keep' in remaining
                count_keep = len([v for v in variant_keep_list if v in variant_remained_list])
                count_delete = len([v for v in variant_delete_list if v in variant_remained_list])
                #
                out_lfn_per_sum_variant_row_dic = {"lfn_per_sum_variant_threshold": lfn_per_sum_variant_threshold,
                           "lfn_read_count_threshold": lfn_read_count_threshold,
                           "count_keep": count_keep, "count_delete": count_delete}
                out_lfn_per_sum_variant_list.append(out_lfn_per_sum_variant_row_dic)
                #
                print(
                    "lfn_per_sum_variant_threshold = {}, lfn_read_count_threshold = {}. count_keep: {}, count_delete: {}"
                    .format(lfn_per_sum_variant_threshold, lfn_read_count_threshold, count_keep, count_delete))
                del (lfn_filter_runner)
                #
                # Increase
                lfn_per_sum_variant_threshold = lfn_per_sum_variant_threshold * 10
                lfn_read_count_threshold = lfn_read_count_threshold * 10
        #
        #
        # lfn_per_sum_variant_threshold = 0.001
        # lfn_read_count_threshold = 10
        # #
        # lfn_filter_runner = FilterLFNRunner(variant_read_count_df)
        # lfn_filter_runner.f7_lfn_delete_absolute_read_count(lfn_read_count_threshold)
        # variant_read_cound_remained_df = lfn_filter_runner.delete_variant_df
        # variant_read_cound_remained_df = variant_read_cound_remained_df.loc[variant_read_cound_remained_df.filter_delete == 0]
        # #
        # variant_remained_list = variant_read_cound_remained_df.variant_id.unique().tolist()
        # # Count how many from 'keep' in remaining
        # count_keep = len([v for v in variant_keep_list if v in variant_remained_list])
        # count_delete = len([v for v in variant_delete_list if v in variant_remained_list])
        # print("lfn_per_sum_variant_threshold = {}, lfn_read_count_threshold = {}. count_keep: {}, count_delete: {}"
        #       .format(lfn_per_sum_variant_threshold, lfn_read_count_threshold, count_keep, count_delete))
        # del(lfn_filter_runner)
        # #
        # #
        # #
        # lfn_per_sum_variant_threshold = 0.01
        # lfn_read_count_threshold = 100
        # #
        # lfn_filter_runner = FilterLFNRunner(variant_read_count_df)
        # lfn_filter_runner.f7_lfn_delete_absolute_read_count(lfn_read_count_threshold)
        # variant_read_cound_remained_df = lfn_filter_runner.delete_variant_df
        # variant_read_cound_remained_df = variant_read_cound_remained_df.loc[variant_read_cound_remained_df.filter_delete == 0]
        # #
        # variant_remained_list = variant_read_cound_remained_df.variant_id.unique().tolist()
        # # Count how many from 'keep' in remaining
        # count_keep = len([v for v in variant_keep_list if v in variant_remained_list])
        # count_delete = len([v for v in variant_delete_list if v in variant_remained_list])
        # print("lfn_per_sum_variant_threshold = {}, lfn_read_count_threshold = {}. count_keep: {}, count_delete: {}"
        #       .format(lfn_per_sum_variant_threshold, lfn_read_count_threshold, count_keep, count_delete))
        # del(lfn_filter_runner)
        # #
        # #
        # #
        # lfn_per_sum_variant_threshold = 1
        # lfn_read_count_threshold = 10000
        # #
        # lfn_filter_runner = FilterLFNRunner(variant_read_count_df)
        # lfn_filter_runner.f7_lfn_delete_absolute_read_count(lfn_read_count_threshold)
        # variant_read_cound_remained_df = lfn_filter_runner.delete_variant_df
        # variant_read_cound_remained_df = variant_read_cound_remained_df.loc[variant_read_cound_remained_df.filter_delete == 0]
        # #
        # variant_remained_list = variant_read_cound_remained_df.variant_id.unique().tolist()
        # # Count how many from 'keep' in remaining
        # count_keep = len([v for v in variant_keep_list if v in variant_remained_list])
        # count_delete = len([v for v in variant_delete_list if v in variant_remained_list])
        # print("lfn_per_sum_variant_threshold = {}, lfn_read_count_threshold = {}. count_keep: {}, count_delete: {}"
        #       .format(lfn_per_sum_variant_threshold, lfn_read_count_threshold, count_keep, count_delete))
        # del(lfn_filter_runner)

