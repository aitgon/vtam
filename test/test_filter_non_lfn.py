import os
import pandas
from unittest import TestCase

from wopmetabarcoding.utils.PathFinder import PathFinder
from wopmetabarcoding.wrapper.FilterLFN import FilterLFNRunner
from wopmetabarcoding.wrapper.FilterNonLFN import FilterNonLFNRunner


class TestFilterNonLFN(TestCase):

    def setUp(self):
        self.__testdir_path = os.path.join(PathFinder.get_module_test_path(), 'test_files')
        #
        self.variant_df = pandas.DataFrame({
            'id':[1,22],
            'sequence_':["tata", "tgtg"],
        })
        #
        # self.MFZR_prerun_COI_corr_all_filters_6
        fin = os.path.join(self.__testdir_path, '6_MFZR_prerun_COI_corr_all_filters.csv')
        df_raw = pandas.read_csv(fin, sep=';', header=0, names=['Variant', 'Biosample', 'Replicate', 'Read_count'])
        df_raw = df_raw[['Variant', 'Biosample', 'Replicate']]
        df_raw[['Biosample2', 'Replicate']] = (df_raw["Replicate"].str.split("-R", n=1, expand=True))
        df_raw['Variant'] = df_raw['Variant'].astype('category')
        df_raw['Biosample'] = df_raw['Biosample'].astype('category')
        df_raw = df_raw[['Variant', 'Biosample', 'Replicate']]
        df_raw['Replicate'] = df_raw['Replicate'].astype('int')
        self.MFZR_prerun_COI_corr_all_filters_6 = pandas.DataFrame({'variant_id': df_raw['Variant'].cat.codes, 'biosample_id': df_raw['Biosample'].cat.codes,
                               'replicate_id': df_raw['Replicate']})
        #
        self.marker_id = 1

    def test_01_f9_delete_min_repln(self):
        #
        filter_non_lfn_runner = FilterNonLFNRunner(self.variant_df, self.MFZR_prerun_COI_corr_all_filters_6, self.marker_id)
        min_repln = 2
        filter_non_lfn_runner.f9_delete_min_repln(min_repln=min_repln)


    # def test_02_f2_f4_lfn_delete_per_sum_variant(self):
    #     lfn_var_threshold = 0.001
    #     self.filter_lfn_runner.f2_f4_lfn_delete_per_sum_variant(lfn_var_threshold)
    #     #
    #     self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 22)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.replicate_id == 1)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.filter_name == 'f2_lfn_delete_per_sum_variant'),
    #                                                                     'filter_delete'].values[0])
    #     self.assertTrue(not self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 22)
    #                                                                  & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
    #                                                                  & (self.filter_lfn_runner.delete_variant_df.replicate_id == 2)
    #                                                                  & (self.filter_lfn_runner.delete_variant_df.filter_name == 'f2_lfn_delete_per_sum_variant'),
    #                                                                     'filter_delete'].values[0])
    #     self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 22)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.replicate_id == 3)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.filter_name == 'f2_lfn_delete_per_sum_variant'),
    #                                                                     'filter_delete'].values[0])
    #
    # def test_03_f2_f4_lfn_delete_per_sum_variant_threshold_specific(self):
    #     lfn_var_threshold = 0.001
    #     lfn_var_threshold_specific = {9: 0.05, 22: 0.01}
    #     self.filter_lfn_runner.f2_f4_lfn_delete_per_sum_variant(lfn_var_threshold, lfn_var_threshold_specific=lfn_var_threshold_specific)
    #     #import pdb; pdb.set_trace()
    #     #
    #     self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 9)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.replicate_id == 1)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.filter_name == 'f5_lfn_var_dep'),
    #                                                                     'filter_delete'].values[0])
    #     self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 9)
    #                                                                  & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
    #                                                                  & (self.filter_lfn_runner.delete_variant_df.replicate_id == 2)
    #                                                                  & (self.filter_lfn_runner.delete_variant_df.filter_name == 'f5_lfn_var_dep'),
    #                                                                     'filter_delete'].values[0])
    #     self.assertTrue(not self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 9)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.biosample_id == 2)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.replicate_id == 1)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.filter_name == 'f5_lfn_var_dep'),
    #                                                                     'filter_delete'].values[0])
    #
    #
    # def test_04_f3_f5_lfn_delete_per_sum_variant_replicate(self):
    #     lfn_var_threshold = 0.005
    #     self.filter_lfn_runner.f3_f5_lfn_delete_per_sum_variant_replicate(lfn_var_threshold)
    #     #
    #     self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[
    #                         (self.filter_lfn_runner.delete_variant_df.variant_id == 12)
    #                         & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
    #                         & (self.filter_lfn_runner.delete_variant_df.replicate_id == 3)
    #                         & (self.filter_lfn_runner.delete_variant_df.filter_name == 'f3_f5_lfn_delete_per_sum_variant_replicate'),
    #                         'filter_delete'].values[0])
    #     self.assertTrue(not self.filter_lfn_runner.delete_variant_df.loc[
    #                         (self.filter_lfn_runner.delete_variant_df.variant_id == 12)
    #                         & (self.filter_lfn_runner.delete_variant_df.biosample_id == 2)
    #                         & (self.filter_lfn_runner.delete_variant_df.replicate_id == 3)
    #                         & (self.filter_lfn_runner.delete_variant_df.filter_name == 'f3_f5_lfn_delete_per_sum_variant_replicate'),
    #                         'filter_delete'].values[0])
    #
    #
    #
    # def test_05_f3_f5_lfn_delete_per_sum_variant_replicate_threshold_specific(self):
    #     lfn_var_threshold = 0.0005
    #     lfn_per_replicate_series_threshold_specific = {9: 0.02, 22: 0.005}
    #     self.filter_lfn_runner.f3_f5_lfn_delete_per_sum_variant_replicate(lfn_var_threshold, lfn_per_replicate_series_threshold_specific=lfn_per_replicate_series_threshold_specific)
    #     #import pdb; pdb.set_trace()
    #     #
    #     self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 22)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.replicate_id == 1)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.filter_name == 'lfn_delete_per_sum_variant_replicate_variant_specific'),
    #                                                                     'filter_delete'].values[0])
    #     self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 9)
    #                                                                  & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
    #                                                                  & (self.filter_lfn_runner.delete_variant_df.replicate_id == 3)
    #                                                                  & (self.filter_lfn_runner.delete_variant_df.filter_name == 'lfn_delete_per_sum_variant_replicate_variant_specific'),
    #                                                                     'filter_delete'].values[0])
    #     self.assertTrue(not self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 9)
    #                                                             & (self.filter_lfn_runner.delete_variant_df.biosample_id == 2)
    #                                                             & (self.filter_lfn_runner.delete_variant_df.replicate_id == 3)
    #                                                             & (self.filter_lfn_runner.delete_variant_df.filter_name == 'lfn_delete_per_sum_variant_replicate_variant_specific'),
    #                                                             'filter_delete'].values[0])
    #     self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 22)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.replicate_id == 3)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.filter_name == 'lfn_delete_per_sum_variant_replicate_variant_specific'),
    #                                                                     'filter_delete'].values[0])
    #
    #
    # def test_06_f7_lfn_delete_absolute_read_count(self):
    #     lfn_read_count_threshold = 10
    #     self.filter_lfn_runner.f7_lfn_delete_absolute_read_count(lfn_read_count_threshold)
    #     #
    #     self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[
    #                         (self.filter_lfn_runner.delete_variant_df.variant_id == 12)
    #                         & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
    #                         & (self.filter_lfn_runner.delete_variant_df.replicate_id == 1)
    #                         & (self.filter_lfn_runner.delete_variant_df.filter_name == 'f7_lfn_delete_absolute_read_count'),
    #                         'filter_delete'].values[0])
    #     self.assertTrue(not self.filter_lfn_runner.delete_variant_df.loc[
    #                         (self.filter_lfn_runner.delete_variant_df.variant_id == 12)
    #                         & (self.filter_lfn_runner.delete_variant_df.biosample_id == 2)
    #                         & (self.filter_lfn_runner.delete_variant_df.replicate_id == 3)
    #                         & (self.filter_lfn_runner.delete_variant_df.filter_name == 'f7_lfn_delete_absolute_read_count'),
    #                         'filter_delete'].values[0])
    #     self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[
    #                         (self.filter_lfn_runner.delete_variant_df.variant_id == 1)
    #                         & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
    #                         & (self.filter_lfn_runner.delete_variant_df.replicate_id == 2)
    #                         & (self.filter_lfn_runner.delete_variant_df.filter_name == 'f7_lfn_delete_absolute_read_count'),
    #                         'filter_delete'].values[0])
    #
    #
    #
    # def test_07_f6_lfn_delete_per_sum_biosample_replicate(self):
    #     lfn_per_replicate_threshold = 0.001
    #
    #     self.filter_lfn_runner.f6_lfn_delete_per_sum_biosample_replicate(lfn_per_replicate_threshold)
    #     #import pdb; pdb.set_trace()
    #     #
    #     self.assertTrue(not self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 9)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.biosample_id == 2)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.replicate_id == 3)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.filter_name == 'f6_lfn_delete_per_sum_biosample_replicate'),
    #                                                                     'filter_delete'].values[0])
    #     self.assertTrue(self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 12)
    #                                                      & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
    #                                                      & (self.filter_lfn_runner.delete_variant_df.replicate_id == 1)
    #                                                      & ( self.filter_lfn_runner.delete_variant_df.filter_name == 'f6_lfn_delete_per_sum_biosample_replicate'),
    #                                                      'filter_delete'].values[0])
    #     self.assertTrue(not self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 12)
    #                                                                  & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
    #                                                                  & (self.filter_lfn_runner.delete_variant_df.replicate_id == 3)
    #                                                                  & (self.filter_lfn_runner.delete_variant_df.filter_name == 'f6_lfn_delete_per_sum_biosample_replicate'),
    #                                                                     'filter_delete'].values[0])
    #
    #     # this cas is false here because we divid 1 by only 46 wich give 0.02>0.001 but in the excel file we divid 1 by  1186
    #     self.assertTrue(not self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 24)
    #                                                             & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
    #                                                             & (self.filter_lfn_runner.delete_variant_df.replicate_id == 3)
    #                                                             & (self.filter_lfn_runner.delete_variant_df.filter_name == 'f6_lfn_delete_per_sum_biosample_replicate'),
    #                                                             'filter_delete'].values[0])
    #       #this cas is false here because we divid 1 by only 161 wich give 0.006>0.001 but in the excel file we divid 1 by  1894
    #     self.assertTrue(not self.filter_lfn_runner.delete_variant_df.loc[(self.filter_lfn_runner.delete_variant_df.variant_id == 24)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.replicate_id == 1)
    #                                                              & (self.filter_lfn_runner.delete_variant_df.filter_name == 'f6_lfn_delete_per_sum_biosample_replicate'),
    #                                                                     'filter_delete'].values[0])

    # def test_05_f5_lfn2_var_dep_mekdad(self):
    #     lfn_per_var = {9: 0.05, 22: 0.01, 0: 0}
    #
    #
    #     self.filter_lfn_runner.f5_lfn2_var_dep_mekdad(self)(lfn_per_var)
    #
    #     self.assertTrue(not self.filter_lfn_runner.delete_variant_df.loc[
    #                         (self.filter_lfn_runner.delete_variant_df.variant_id == 9)
    #                         & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
    #                         & (self.filter_lfn_runner.delete_variant_df.replicate_id == 2)
    #                         & (self.filter_lfn_runner.delete_variant_df.filter_name == 'f5_lfn2_var_dep_mekdad'),
    #                         'filter_passed'].values[0])
    #     self.assertTrue(not self.filter_lfn_runner.delete_variant_df.loc[
    #                         (self.filter_lfn_runner.delete_variant_df.variant_id == 22)
    #                         & (self.filter_lfn_runner.delete_variant_df.biosample_id == 1)
    #                         & (self.filter_lfn_runner.delete_variant_df.replicate_id == 2)
    #                         & (self.filter_lfn_runner.delete_variant_df.filter_name == 'f5_lfn2_var_dep_mekdad'),
    #                         'filter_passed'].values[0] == False)