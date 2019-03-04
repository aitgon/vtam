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
        fin1 = os.path.join(self.__testdir_path, '7_MFZR_prerun_COI_corr_WO_incoherent_var.csv')
        df_raw = pandas.read_csv(fin, sep=';', header=0, names=['Variant', 'Biosample', 'Replicate', 'Read_count'])
        df_raw1 = pandas.read_csv(fin1, sep=';', header=0, names=['Variant', 'Biosample', 'Replicate', 'Read_count'])
        df_raw = df_raw[['Variant', 'Biosample', 'Replicate']]
        df_raw1 = df_raw1[['Variant', 'Biosample', 'Replicate']]
        # df_raw[['Biosample2', 'Replicate']] = (df_raw["Replicate"].str.split("-R", n=1, expand=True))
        df_raw['Variant'] = df_raw['Variant'].astype('category')
        df_raw['Biosample'] = df_raw['Biosample'].astype('category')
        # import pdb; pdb.set_trace()
        df_raw = df_raw[['Variant', 'Biosample', 'Replicate']]
        # df_raw['Replicate'] = df_raw['Replicate'].astype('int')
        self.MFZR_prerun_COI_corr_all_filters_6 = pandas.DataFrame({'variant_id': df_raw['Variant'].cat.codes, 'biosample_id': df_raw['Biosample'].cat.codes,
                               'replicate_id': df_raw['Replicate']})
        #
        self.marker_id = 1

    def test_01_f9_delete_min_repln(self):
        #
        filter_non_lfn_runner = FilterNonLFNRunner(self.variant_df, self.MFZR_prerun_COI_corr_all_filters_6, self.marker_id)
        min_repln = 2
        var_prop = 0
        pcr_error_by_sample = True
        # f9_pcr_error(var_prop=var_prop, pcr_error_by_sample=pcr_error_by_sample)
        filter_non_lfn_runner.f9_delete_min_replicate_number(min_replicate_number=min_repln)
        #
        filter_non_lfn_runner.delete_variant_df.loc[filter_non_lfn_runner.delete_variant_df.filter_delete == 0].shape
        #
        # Nb of rows in wopmetabarcodin/test/test_files/7_MFZR_prerun_COI_corr_WO_incoherent_var.csv
        # 836
        nb_variant_biosample_replicates_that_passed_the_filter\
            = filter_non_lfn_runner.delete_variant_df.loc[filter_non_lfn_runner.delete_variant_df.filter_delete==0].shape[0]
        self.assertTrue(nb_variant_biosample_replicates_that_passed_the_filter==836)
