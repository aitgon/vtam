import os

import pandas
from unittest import TestCase

from wopmetabarcoding.utils.PathManager import PathFinder
from wopmetabarcoding.wrapper.FilterLFNutilities import FilterLFNRunner, f1_lfn_delete_singleton


class TestSingleton(TestCase):
    def setUp(self):
        self.variant_df = pandas.DataFrame({
            'id': [1, 22],
            'sequence_': ["tata", "tgtg"],
        })
        self.variant_read_count_df = pandas.DataFrame({
            'run_id': 150*[1],
            'marker_id': 150*[1],
            'variant_sequence': [1]*6 + [2]*6 + [3]*6 + [4]*6 + [5]*6 + [6]*6 + [7]*6 + [8]*6 + [9]*6
                                + [10]*6 + [11]*6 + [12]*6 + [13]*6 +[14]*6 + [15]*6 + [16]*6 + [17]*6 + [18]*6 + [19]*6
                                + [20]*6 + [21]*6 + [22]*6 + [23]*6 + [24]*6 + [25]*6,
            'biosample_id': 25*(3*[1] + 3*[2]),
            'replicate_id': 50* [1,2,3],
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


    def test01_lfn_delete_singleton(self):
        passed_lfn_delete_singleton_df = f1_lfn_delete_singleton(self.variant_read_count_df)
        #
        nb_variant_that_passed_the_filter= passed_lfn_delete_singleton_df.shape[0]
        self.assertTrue(nb_variant_that_passed_the_filter == 126)


class TestFilterLFN(TestCase):

    def setUp(self):
        #
        self.__testdir_path = os.path.join(PathFinder.get_module_test_path())
        self.lfn_var_threshold_specific = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "test_files", "lfn_var_threshold_specific.tsv")

        #
        self.variant_df = pandas.DataFrame({
            'id':[1,22],
            'sequence_':["tata", "tgtg"],
        })
        self.variant_read_count_df = pandas.DataFrame({
            'run_id': [1]*150,
            'marker_id': 150*[1],
            'variant_id': [1]*6 + [2]*6 + [3]*6 + [4]*6 + [5]*6 + [6]*6+ [7]*6 + [8]*6 + [9]*6 + [10]*6 + [11]*6 + [12]*6 +  [13]*6+
                          [14]*6 + [15]*6 + [16]*6 + [17]*6 + [18]*6 + [19]*6 + [20]*6+ [21]*6 + [22]*6 + [23]*6 + [24]*6+ [25]*6,
            'biosample_id':[1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,
                            1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,
                            1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2],
            'replicate_id':[1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,
                            1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,
                            1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3],
            'read_count':[
                10,5,0,249,58,185,
                68,54,100,0,0,0,
                0,0,0,258,126,500,
                0,0,0,0,1,0,
                0,0,1,0,0,0,
                1524,1815,789,118,98,50,
                1,0,0,0,0,0,
                0,1,0,0,0,0,
                125,214,20,1284,1789,1913,
                0,1,0,0,1,0,
                15,0,1,0,0,25,
                0,0,2,598,50,875,
                2,60,12,1,0,0,
                1,0,0,0,0,2,
                0,3,0,0,5,0,
                65,98,152,2,0,1,
                52,74,85,0,0,0,
                1,0,0,5,0,8,
                5,0,1,0,0,21,
                0,0,0,524,658,125,
                0,0,0,2,0,10,
                25,58,23,10980,8999,13814,
                0,5,0,0,2,0,
                1,0,1,1,0,284,
                0,2,0,0,5,0,
                  ],
        })
        self.marker_id = 1
        #
        self.filter_lfn_runner = FilterLFNRunner(self.variant_read_count_df)

    def test_02_f2_f4_lfn_delete_variant(self):
        lfn_per_variant_threshold = 0.001
        self.filter_lfn_runner.f2_f4_lfn_delete_variant(lfn_variant_threshold=lfn_per_variant_threshold)
        #
        self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 22)
                                                                     & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
                                                                     & (self.filter_lfn_runner.delete_variant_df.replicate_id == 1)
                                                                     & (self.filter_lfn_runner.delete_variant_df.filter_id == 2),
                                                                        'filter_delete'].values[0])
        self.assertTrue(not self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 22)
                                                                         & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
                                                                         & (self.filter_lfn_runner.delete_variant_df.replicate_id == 2)
                                                                         & (self.filter_lfn_runner.delete_variant_df.filter_id == 2),
                                                                        'filter_delete'].values[0])
        self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 22)
                                                                     & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
                                                                     & (self.filter_lfn_runner.delete_variant_df.replicate_id == 3)
                                                                     & (self.filter_lfn_runner.delete_variant_df.filter_id == 2),
                                                                        'filter_delete'].values[0])

    def test_03_f2_f4_lfn_delete_variant_threshold_specific(self):
        lfn_var_threshold = 0.001
        # lfn_var_threshold_specific = {9: 0.05, 22: 0.01}
        # input tsv file of threshold specific and create a dictionnary
        lfn_var_threshold_specific_df = pandas.read_csv(self.lfn_var_threshold_specific, sep='\t', header=0)
        # lfn_var_threshold_specific={}
        # for  row in lfn_var_threshold_specific_df.itertuples():
        #    lfn_var_threshold_specific[row.variant_id]= float(row.variant_id_threshold)
        # import pdb;
        # pdb.set_trace()

        self.filter_lfn_runner.f2_f4_lfn_delete_variant(lfn_var_threshold, threshold_specific_df=lfn_var_threshold_specific_df)
        #
        self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 9)
                                                                     & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
                                                                     & (self.filter_lfn_runner.delete_variant_df.replicate_id == 1)
                                                                     & (self.filter_lfn_runner.delete_variant_df.filter_id == 4),
                                                                        'filter_delete'].values[0])
        self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 9)
                                                                     & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
                                                                     & (self.filter_lfn_runner.delete_variant_df.replicate_id == 2)
                                                                     & (self.filter_lfn_runner.delete_variant_df.filter_id == 4),
                                                                        'filter_delete'].values[0])
        self.assertTrue(not self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 9)
                                                                         & (self.filter_lfn_runner.delete_variant_df.biosample_id == 2)
                                                                         & (self.filter_lfn_runner.delete_variant_df.replicate_id == 1)
                                                                         & (self.filter_lfn_runner.delete_variant_df.filter_id == 4),
                                                                        'filter_delete'].values[0])


    def test_04_f3_f5_lfn_delete_per_sum_variant_replicate(self):
        lfn_per_replicate_threshold = 0.005
        self.filter_lfn_runner.f3_f5_lfn_delete_variant_replicate(
            lfn_variant_replicate_threshold=lfn_per_replicate_threshold)
        #
        self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[
                            (self.filter_lfn_runner.delete_variant_df.variant_id == 12)
                            & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
                            & (self.filter_lfn_runner.delete_variant_df.replicate_id == 3)
                            & (self.filter_lfn_runner.delete_variant_df.filter_id == 3),
                            'filter_delete'].values[0])
        self.assertTrue(not self.filter_lfn_runner.delete_variant_df.loc[
                            (self.filter_lfn_runner.delete_variant_df.variant_id == 12)
                            & (self.filter_lfn_runner.delete_variant_df.biosample_id == 2)
                            & (self.filter_lfn_runner.delete_variant_df.replicate_id == 3)
                            & (self.filter_lfn_runner.delete_variant_df.filter_id == 3),
                            'filter_delete'].values[0])



    def test_05_f3_f5_lfn_delete_per_sum_variant_replicate_threshold_specific(self):
        lfn_var_threshold = 0.0005
        # lfn_per_replicate_series_threshold_specific = {9: 0.02, 22: 0.005}
        lfn_var_threshold_specific_df = pandas.read_csv(self.lfn_var_threshold_specific, sep='\t', header=0)
        self.filter_lfn_runner.f3_f5_lfn_delete_variant_replicate(lfn_var_threshold, threshold_specific_df=lfn_var_threshold_specific_df)
        #import pdb; pdb.set_trace()
        #
        self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 22)
                                                                     & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
                                                                     & (self.filter_lfn_runner.delete_variant_df.replicate_id == 1)
                                                                     & (self.filter_lfn_runner.delete_variant_df.filter_id == 5),
                                                                        'filter_delete'].values[0])
        self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 9)
                                                                     & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
                                                                     & (self.filter_lfn_runner.delete_variant_df.replicate_id == 3)
                                                                     & (self.filter_lfn_runner.delete_variant_df.filter_id == 5),
                                                                        'filter_delete'].values[0])
        self.assertTrue(not self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 9)
                                                                         & (self.filter_lfn_runner.delete_variant_df.biosample_id == 2)
                                                                         & (self.filter_lfn_runner.delete_variant_df.replicate_id == 3)
                                                                         & (self.filter_lfn_runner.delete_variant_df.filter_id == 5),
                                                                'filter_delete'].values[0])
        self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 22)
                                                                     & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
                                                                     & (self.filter_lfn_runner.delete_variant_df.replicate_id == 3)
                                                                     & (self.filter_lfn_runner.delete_variant_df.filter_id == 5),
                                                                        'filter_delete'].values[0])


    def test_06_f7_lfn_delete_absolute_read_count(self):
        lfn_read_count_threshold = 10
        self.filter_lfn_runner.f7_lfn_delete_absolute_read_count(lfn_read_count_threshold)
        #
        self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[
                            (self.filter_lfn_runner.delete_variant_df.variant_id == 12)
                            & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
                            & (self.filter_lfn_runner.delete_variant_df.replicate_id == 1)
                            & (self.filter_lfn_runner.delete_variant_df.filter_id == 7),
                            'filter_delete'].values[0])
        self.assertTrue(not self.filter_lfn_runner.delete_variant_df.loc[
                            (self.filter_lfn_runner.delete_variant_df.variant_id == 12)
                            & (self.filter_lfn_runner.delete_variant_df.biosample_id == 2)
                            & (self.filter_lfn_runner.delete_variant_df.replicate_id == 3)
                            & (self.filter_lfn_runner.delete_variant_df.filter_id == 7),
                            'filter_delete'].values[0])
        self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[
                            (self.filter_lfn_runner.delete_variant_df.variant_id == 1)
                            & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
                            & (self.filter_lfn_runner.delete_variant_df.replicate_id == 2)
                            & (self.filter_lfn_runner.delete_variant_df.filter_id == 7),
                            'filter_delete'].values[0])



    def test_07_f6_lfn_delete_per_sum_biosample_replicate(self):
        lfn_per_replicate_threshold = 0.001

        self.filter_lfn_runner.f6_lfn_delete_biosample_replicate(lfn_per_replicate_threshold)
        #import pdb; pdb.set_trace()
        #
        self.assertTrue(not self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 9)
                                                                         & (self.filter_lfn_runner.delete_variant_df.biosample_id == 2)
                                                                         & (self.filter_lfn_runner.delete_variant_df.replicate_id == 3)
                                                                         & (self.filter_lfn_runner.delete_variant_df.filter_id == 6),
                                                                        'filter_delete'].values[0])
        self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 12)
                                                                     & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
                                                                     & (self.filter_lfn_runner.delete_variant_df.replicate_id == 1)
                                                                     & (self.filter_lfn_runner.delete_variant_df.filter_id == 6),
                                                         'filter_delete'].values[0])
        self.assertTrue(not self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 12)
                                                                         & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
                                                                         & (self.filter_lfn_runner.delete_variant_df.replicate_id == 3)
                                                                         & (self.filter_lfn_runner.delete_variant_df.filter_id == 6),
                                                                        'filter_delete'].values[0])


        self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 24)
                                                                         & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
                                                                         & (self.filter_lfn_runner.delete_variant_df.replicate_id == 3)
                                                                         & (self.filter_lfn_runner.delete_variant_df.filter_id == 6),
                                                                'filter_delete'].values[0])

        self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 24)
                                                                         & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
                                                                         & (self.filter_lfn_runner.delete_variant_df.replicate_id == 1)
                                                                         & (self.filter_lfn_runner.delete_variant_df.filter_id == 6),
                                                                        'filter_delete'].values[0])



    def test_08_f8_lfn_delete_do_not_pass_all_filters(self):
        lfn_per_variant_threshold = 0.001
        self.filter_lfn_runner.f2_f4_lfn_delete_variant(lfn_variant_threshold=lfn_per_variant_threshold)
        lfn_per_replicate_threshold = 0.005
        self.filter_lfn_runner.f3_f5_lfn_delete_variant_replicate(
            lfn_variant_replicate_threshold=lfn_per_replicate_threshold)
        #
        self.filter_lfn_runner.f8_lfn_delete_do_not_pass_all_filters()
        #
        self.assertTrue(not self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.run_id == 1)
                                                                         & (self.filter_lfn_runner.delete_variant_df.variant_id == 1)
                                                                         & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
                                                                         & (self.filter_lfn_runner.delete_variant_df.replicate_id == 1)
                                                                         & (self.filter_lfn_runner.delete_variant_df.filter_id == 8),
                                                                        'filter_delete'].values[0])
        self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.run_id == 1)                                                                         & (self.filter_lfn_runner.delete_variant_df.variant_id == 1)
                                                                         & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
                                                                         & (self.filter_lfn_runner.delete_variant_df.replicate_id == 3)
                                                                         & (self.filter_lfn_runner.delete_variant_df.filter_id == 8),
                                                                        'filter_delete'].values[0])


