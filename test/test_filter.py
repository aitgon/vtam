import pandas
from unittest import TestCase
from wopmetabarcoding.wrapper.LFNFilters import LFNFilterRunner

class TestFilter(TestCase):

    def setUp(self):
        self.variant_df = pandas.DataFrame({
            'id':[1,22],
            'sequence_':["tata", "tgtg"],
        })
        self.variant_read_count_df = pandas.DataFrame({
            'variant_id': [1]*6 + [12]*6 + [9]*6 +[22]*6 +[24]*6,
            'biosample_id':[1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2],
            'replicate_id':[1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3],
            'read_count':[
                10,5,0,249,58,185,
                0,0,2,598,50,875,
                125,214,20,1284,1798,1913,
                25,58,23,10980,8999,13814,
                1,0,1,1,0,284,
            ],
        })
        self.marker_id = 1
        #
        self.filter_runner = LFNFilterRunner(self.variant_df, self.variant_read_count_df, self.marker_id)

    def test_02_f2_f4_lfn_delete_per_sum_variant(self):
        lfn_var_threshold = 0.001
        self.filter_runner.f2_f4_lfn_delete_per_sum_variant(lfn_var_threshold)
        #
        self.assertTrue(self.filter_runner.delete_variant_df.loc[(self.filter_runner.delete_variant_df.variant_id == 22)
                                                                 & (self.filter_runner.delete_variant_df.biosample_id == 1)
                                                                 & (self.filter_runner.delete_variant_df.replicate_id == 1)
                                                                 & (self.filter_runner.delete_variant_df.filter_name == 'f2_lfn_var'),
                                                                        'filter_delete'].values[0])
        self.assertTrue(not self.filter_runner.delete_variant_df.loc[(self.filter_runner.delete_variant_df.variant_id == 22)
                                                                     & (self.filter_runner.delete_variant_df.biosample_id == 1)
                                                                     & (self.filter_runner.delete_variant_df.replicate_id == 2)
                                                                     & (self.filter_runner.delete_variant_df.filter_name == 'f2_lfn_var'),
                                                                        'filter_delete'].values[0])
        self.assertTrue(self.filter_runner.delete_variant_df.loc[(self.filter_runner.delete_variant_df.variant_id == 22)
                                                                 & (self.filter_runner.delete_variant_df.biosample_id == 1)
                                                                 & (self.filter_runner.delete_variant_df.replicate_id == 3)
                                                                 & (self.filter_runner.delete_variant_df.filter_name == 'f2_lfn_var'),
                                                                        'filter_delete'].values[0])

    def test_03_f2_f4_lfn_delete_per_sum_variant_threshold_specific(self):
        lfn_var_threshold = 0.001
        lfn_var_threshold_specific = {9: 0.05, 22: 0.01}
        self.filter_runner.f2_f4_lfn_delete_per_sum_variant(lfn_var_threshold, lfn_var_threshold_specific=lfn_var_threshold_specific)
        #import pdb; pdb.set_trace()
        #
        self.assertTrue(self.filter_runner.delete_variant_df.loc[(self.filter_runner.delete_variant_df.variant_id == 9)
                                                                 & (self.filter_runner.delete_variant_df.biosample_id == 1)
                                                                 & (self.filter_runner.delete_variant_df.replicate_id == 1)
                                                                 & (self.filter_runner.delete_variant_df.filter_name == 'f5_lfn_var_dep'),
                                                                        'filter_delete'].values[0])
        self.assertTrue(self.filter_runner.delete_variant_df.loc[(self.filter_runner.delete_variant_df.variant_id == 9)
                                                                     & (self.filter_runner.delete_variant_df.biosample_id == 1)
                                                                     & (self.filter_runner.delete_variant_df.replicate_id == 2)
                                                                     & (self.filter_runner.delete_variant_df.filter_name == 'f5_lfn_var_dep'),
                                                                        'filter_delete'].values[0])
        self.assertTrue(not self.filter_runner.delete_variant_df.loc[(self.filter_runner.delete_variant_df.variant_id == 9)
                                                                 & (self.filter_runner.delete_variant_df.biosample_id == 2)
                                                                 & (self.filter_runner.delete_variant_df.replicate_id == 1)
                                                                 & (self.filter_runner.delete_variant_df.filter_name == 'f5_lfn_var_dep'),
                                                                        'filter_delete'].values[0])


    def test_04_f3_f5_lfn_delete_per_sum_variant_replicate(self):
        lfn_var_threshold = 0.005
        self.filter_runner.f3_f5_lfn_delete_per_sum_variant_replicate(lfn_var_threshold)
        #
        self.assertTrue(self.filter_runner.delete_variant_df.loc[
                            (self.filter_runner.delete_variant_df.variant_id == 12)
                            & (self.filter_runner.delete_variant_df.biosample_id == 1)
                            & (self.filter_runner.delete_variant_df.replicate_id == 3)
                            & (self.filter_runner.delete_variant_df.filter_name == 'f3_f5_lfn_delete_per_sum_variant_replicate'),
                            'filter_delete'].values[0])
        self.assertTrue(not self.filter_runner.delete_variant_df.loc[
                            (self.filter_runner.delete_variant_df.variant_id == 12)
                            & (self.filter_runner.delete_variant_df.biosample_id == 2)
                            & (self.filter_runner.delete_variant_df.replicate_id == 3)
                            & (self.filter_runner.delete_variant_df.filter_name == 'f3_f5_lfn_delete_per_sum_variant_replicate'),
                            'filter_delete'].values[0])



    def test_05_f3_f5_lfn_delete_per_sum_variant_replicate_threshold_specific(self):
        lfn_var_threshold = 0.0005
        lfn_per_replicate_series_threshold_specific = {9: 0.02, 22: 0.005}
        self.filter_runner.f3_f5_lfn_delete_per_sum_variant_replicate(lfn_var_threshold, lfn_per_replicate_series_threshold_specific=lfn_per_replicate_series_threshold_specific)
        #import pdb; pdb.set_trace()
        #
        self.assertTrue(self.filter_runner.delete_variant_df.loc[(self.filter_runner.delete_variant_df.variant_id == 22)
                                                                 & (self.filter_runner.delete_variant_df.biosample_id == 1)
                                                                 & (self.filter_runner.delete_variant_df.replicate_id == 1)
                                                                 & (self.filter_runner.delete_variant_df.filter_name == 'lfn_delete_per_sum_variant_replicate_variant_specific'),
                                                                        'filter_delete'].values[0])
        self.assertTrue(self.filter_runner.delete_variant_df.loc[(self.filter_runner.delete_variant_df.variant_id == 9)
                                                                     & (self.filter_runner.delete_variant_df.biosample_id == 1)
                                                                     & (self.filter_runner.delete_variant_df.replicate_id == 3)
                                                                     & (self.filter_runner.delete_variant_df.filter_name == 'lfn_delete_per_sum_variant_replicate_variant_specific'),
                                                                        'filter_delete'].values[0])
        self.assertTrue(not self.filter_runner.delete_variant_df.loc[(self.filter_runner.delete_variant_df.variant_id == 9)
                                                                & (self.filter_runner.delete_variant_df.biosample_id == 2)
                                                                & (self.filter_runner.delete_variant_df.replicate_id == 3)
                                                                & (self.filter_runner.delete_variant_df.filter_name == 'lfn_delete_per_sum_variant_replicate_variant_specific'),
                                                                'filter_delete'].values[0])
        self.assertTrue(self.filter_runner.delete_variant_df.loc[(self.filter_runner.delete_variant_df.variant_id == 22)
                                                                 & (self.filter_runner.delete_variant_df.biosample_id == 1)
                                                                 & (self.filter_runner.delete_variant_df.replicate_id == 3)
                                                                 & (self.filter_runner.delete_variant_df.filter_name == 'lfn_delete_per_sum_variant_replicate_variant_specific'),
                                                                        'filter_delete'].values[0])


    def test_06_f7_lfn_delete_absolute_read_count(self):
        lfn_read_count_threshold = 10
        self.filter_runner.f7_lfn_delete_absolute_read_count(lfn_read_count_threshold)
        #
        self.assertTrue(self.filter_runner.delete_variant_df.loc[
                            (self.filter_runner.delete_variant_df.variant_id == 12)
                            & (self.filter_runner.delete_variant_df.biosample_id == 1)
                            & (self.filter_runner.delete_variant_df.replicate_id == 1)
                            & (self.filter_runner.delete_variant_df.filter_name == 'f7_lfn_delete_absolute_read_count'),
                            'filter_delete'].values[0])
        self.assertTrue(not self.filter_runner.delete_variant_df.loc[
                            (self.filter_runner.delete_variant_df.variant_id == 12)
                            & (self.filter_runner.delete_variant_df.biosample_id == 2)
                            & (self.filter_runner.delete_variant_df.replicate_id == 3)
                            & (self.filter_runner.delete_variant_df.filter_name == 'f7_lfn_delete_absolute_read_count'),
                            'filter_delete'].values[0])
        self.assertTrue(self.filter_runner.delete_variant_df.loc[
                            (self.filter_runner.delete_variant_df.variant_id == 1)
                            & (self.filter_runner.delete_variant_df.biosample_id == 1)
                            & (self.filter_runner.delete_variant_df.replicate_id == 2)
                            & (self.filter_runner.delete_variant_df.filter_name == 'f7_lfn_delete_absolute_read_count'),
                            'filter_delete'].values[0])



    # def test_07_f1_lfn1_per_replicate_new(self):
    #     lfn_per_replicate_threshold = 0.001
    #
    #     self.filter_runner.f6_lfn_delete_per_sum_biosample_replicate(lfn_per_replicate_threshold)
    #     #import pdb; pdb.set_trace()
    #     #
    #     self.assertTrue(not self.filter_runner.delete_variant_df.loc[(self.filter_runner.delete_variant_df.variant_id == 9)
    #                                                              & (self.filter_runner.delete_variant_df.biosample_id == 2)
    #                                                              & (self.filter_runner.delete_variant_df.replicate_id == 3)
    #                                                              & (self.filter_runner.delete_variant_df.filter_name == 'f1_lfn1_per_replicate_new_07'),
    #                                                                     'filter_delete'].values[0])
    #     self.assertTrue(self.filter_runner.delete_variant_df.loc[(self.filter_runner.delete_variant_df.variant_id == 12)
    #                                                      & (self.filter_runner.delete_variant_df.biosample_id == 1)
    #                                                      & (self.filter_runner.delete_variant_df.replicate_id == 1)
    #                                                      & ( self.filter_runner.delete_variant_df.filter_name == 'f1_lfn1_per_replicate_new_07'),
    #                                                      'filter_delete'].values[0])
    #     self.assertTrue(not self.filter_runner.delete_variant_df.loc[(self.filter_runner.delete_variant_df.variant_id == 12)
    #                                                                  & (self.filter_runner.delete_variant_df.biosample_id == 1)
    #                                                                  & (self.filter_runner.delete_variant_df.replicate_id == 3)
    #                                                                  & (self.filter_runner.delete_variant_df.filter_name == 'f1_lfn1_per_replicate_new_07'),
    #                                                                     'filter_delete'].values[0])
    #
    #     # this cas is false here because we divid 1 by only 46 wich give 0.02>0.001 but in the excel file we divid 1 by  1186
    #     self.assertTrue(not self.filter_runner.delete_variant_df.loc[(self.filter_runner.delete_variant_df.variant_id == 24)
    #                                                             & (self.filter_runner.delete_variant_df.biosample_id == 1)
    #                                                             & (self.filter_runner.delete_variant_df.replicate_id == 3)
    #                                                             & (self.filter_runner.delete_variant_df.filter_name == 'f1_lfn1_per_replicate_new_07'),
    #                                                             'filter_delete'].values[0])
    #       #this cas is false here because we divid 1 by only 161 wich give 0.006>0.001 but in the excel file we divid 1 by  1894
    #     self.assertTrue(not self.filter_runner.delete_variant_df.loc[(self.filter_runner.delete_variant_df.variant_id == 24)
    #                                                              & (self.filter_runner.delete_variant_df.biosample_id == 1)
    #                                                              & (self.filter_runner.delete_variant_df.replicate_id == 1)
    #                                                              & (self.filter_runner.delete_variant_df.filter_name == 'f1_lfn1_per_replicate_new_07'),
    #                                                                     'filter_delete'].values[0])

    # def test_05_f5_lfn2_var_dep_mekdad(self):
    #     lfn_per_var = {9: 0.05, 22: 0.01, 0: 0}
    #
    #
    #     self.filter_runner.f5_lfn2_var_dep_mekdad(self)(lfn_per_var)
    #
    #     self.assertTrue(not self.filter_runner.delete_variant_df.loc[
    #                         (self.filter_runner.delete_variant_df.variant_id == 9)
    #                         & (self.filter_runner.delete_variant_df.biosample_id == 1)
    #                         & (self.filter_runner.delete_variant_df.replicate_id == 2)
    #                         & (self.filter_runner.delete_variant_df.filter_name == 'f5_lfn2_var_dep_mekdad'),
    #                         'filter_passed'].values[0])
    #     self.assertTrue(not self.filter_runner.delete_variant_df.loc[
    #                         (self.filter_runner.delete_variant_df.variant_id == 22)
    #                         & (self.filter_runner.delete_variant_df.biosample_id == 1)
    #                         & (self.filter_runner.delete_variant_df.replicate_id == 2)
    #                         & (self.filter_runner.delete_variant_df.filter_name == 'f5_lfn2_var_dep_mekdad'),
    #                         'filter_passed'].values[0] == False)