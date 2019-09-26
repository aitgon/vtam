import os
import sqlite3
from unittest import TestCase

import pandas

from wopmetabarcoding.utils.PathFinder import PathFinder
from wopmetabarcoding.utils.utilities import create_step_tmp_dir
from wopmetabarcoding.wrapper.FilterLFNutilities import FilterLFNRunner
from wopmetabarcoding.wrapper.FilterMinReplicateNumber import f9_delete_min_replicate_number


class TestOptimizeF7(TestCase):

    def setUp(self):
        self.__testdir_path = os.path.join(PathFinder.get_module_test_path())
        self.variant_read_count_path = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "test_files",
                                                  "optimize_f7", "variant_read_count.tsv")
        self.variants_optimize_path = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "test_files",
                                                  "optimize_f7", "variants_optimize.tsv")
        self.this_step_tmp_dir = create_step_tmp_dir(__file__)
        #

    def test_01(self):


        ####
        #
        #  read_count of keep and delete
        #
        ####
        variant_read_count_df = pandas.read_csv(self.variant_read_count_path, sep='\t', header=0)
        #
        variants_optimize_df = pandas.read_csv(self.variants_optimize_path, sep='\t', header=0)
        #
        # variant_keep_list = (
        # variant_read_count_df.loc[variant_read_count_df.action == 'keep']).variant_id.unique().tolist()
        # variant_delete_list = (
        # variant_read_count_df.loc[variant_read_count_df.action == 'delete']).variant_id.unique().tolist()


        variant_keep_df = variants_optimize_df.loc[variants_optimize_df.action == 'keep']
        variant_delete_negative_df = variants_optimize_df.loc[(variants_optimize_df.action == 'delete') &
                                                     (variants_optimize_df.biosample_type == 'negative')]
        variant_delete_real_df = variants_optimize_df.loc[(variants_optimize_df.action == 'delete') &
                                                     (variants_optimize_df.biosample_type == 'real')]

        ##############
        #
        # Main loop through parameter values
        #
        ##############
        #
        out_lfn_variant_list = []
        #
        min_repln = 2
        lfn_biosample_replicate_threshold = 0.001
        #
        lfn_variant_threshold = 0.001 # default value
        lfn_variant_threshold_range = [i / 10000 for i in range(1, 1000)]
        lfn_read_count_threshold = 10
        lfn_read_count_threshold_range = [10*i for i in range(1, 1000)]
        # lfn_variant_threshold_max = lfn_variant_threshold_min * 1000
        # lfn_variant_threshold = lfn_variant_threshold_min
        #
        count_keep = 0
        count_keep_max = 0
        #
        while count_keep >= count_keep_max:
            #
            # lfn_read_count_threshold_min = 10
            # lfn_read_count_threshold_max = lfn_read_count_threshold_min * 1000
            # lfn_read_count_threshold = lfn_read_count_threshold_min
            #
            # while lfn_read_count_threshold <= lfn_read_count_threshold_max and count_keep >= count_keep_max:
                #
            lfn_filter_runner = FilterLFNRunner(variant_read_count_df)

            ###################
            #
            # Filter lfn_variant
            #
            ####################

            lfn_filter_runner.f2_f4_lfn_delete_variant(lfn_variant_threshold)

            ###################
            #
            # Filter lfn_biosample_replicate
            #
            ####################

            lfn_filter_runner.f6_lfn_delete_biosample_replicate(lfn_biosample_replicate_threshold)

            ###################
            #
            # Filter absolute read count
            #
            ####################

            lfn_filter_runner.f7_lfn_delete_absolute_read_count(lfn_read_count_threshold)

            ###################
            #
            # f8_lfn_delete_do_not_pass_all_filters
            #
            ####################

            lfn_filter_runner.f8_lfn_delete_do_not_pass_all_filters()

            variant_read_count_remained_df = lfn_filter_runner.delete_variant_df

            variant_read_count_remained_df = variant_read_count_remained_df.loc[
                (variant_read_count_remained_df.filter_id == 8) &
                (variant_read_count_remained_df.filter_delete == 0)]

            ##########################################################
            #
            # f9_delete_min_replicate_number
            #
            ##########################################################

            variant_read_count_remained_df = f9_delete_min_replicate_number(variant_read_count_remained_df, min_repln)
            variant_read_count_remained_df = variant_read_count_remained_df.loc[
                (variant_read_count_remained_df.filter_delete == 0)]
            variant_read_count_remained_df.drop('filter_delete', axis=1, inplace=True)

            ##########################################################
            #
            # Count keep
            #
            ##########################################################

            variant_read_count_remained_df = variant_read_count_remained_df[['run_id', 'marker_id', 'variant_id', 'biosample_id']]
            variant_read_count_remained_df.drop_duplicates(inplace=True)
            variant_read_count_remained_keep_df = variant_read_count_remained_df.merge(variant_keep_df,
                                                 on=['run_id', 'marker_id', 'variant_id', 'biosample_id']).drop_duplicates()
            count_keep = variant_read_count_remained_keep_df.shape[0]

            ##########################################################
            #
            # Count delete
            #
            ##########################################################

            variant_read_count_remained_delete_negative_df = variant_read_count_remained_df.merge(
                variant_delete_negative_df, on=['run_id', 'marker_id', 'biosample_id']).drop_duplicates()
            variant_read_count_remained_delete_real_df = variant_read_count_remained_df.merge(
                variant_delete_real_df, on=['run_id', 'marker_id', 'variant_id', 'biosample_id']).drop_duplicates()
            count_delete = variant_read_count_remained_delete_negative_df.shape[0] \
                           + variant_read_count_remained_delete_real_df.shape[0]

            # variant_read_count_remained_df[
            #     ['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id', 'filter_delete']].groupby(
            #     ['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id']).sum().reset_index()
            #
            # variant_read_count_remained_df = variant_read_count_remained_df.loc[
            #     variant_read_count_remained_df.filter_delete == 0]
            # #
            # variant_remained_list = variant_read_count_remained_df.variant_id.unique().tolist()
            # #  Count how many from 'keep' in remaining
            # count_keep = len([v for v in variant_keep_list if v in variant_remained_list])
            # count_delete = len([v for v in variant_delete_list if v in variant_remained_list])
            # #
            out_lfn_variant_row_dic = {"lfn_variant_threshold": lfn_variant_threshold,
                       "lfn_read_count_threshold": lfn_read_count_threshold,
                       "count_keep": count_keep, "count_delete": count_delete}
            out_lfn_variant_list.append(out_lfn_variant_row_dic)
            del (lfn_filter_runner)
            #
            print(out_lfn_variant_row_dic, count_keep_max)
            # if count_keep < count_keep_max:
            #     break
            if count_keep > count_keep_max:
                count_keep_max = count_keep
            #
            # Increase
            lfn_read_count_threshold = lfn_read_count_threshold + 5
            lfn_variant_threshold = lfn_variant_threshold + 0.0005
        # import pdb; pdb.set_trace()
        self.assertTrue(out_lfn_variant_list[0]
                        == {'lfn_variant_threshold': 0.001, 'lfn_read_count_threshold': 10, 'count_keep': 12,
                            'count_delete': 3})
        out_lfn_variant_df = pandas.DataFrame(out_lfn_variant_list) # output
        # print(pandas.DataFrame(out_lfn_variant_list))
