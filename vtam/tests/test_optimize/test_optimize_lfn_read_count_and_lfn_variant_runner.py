import os
import pandas
import unittest

from vtam.utils.OptimizeLFNreadCountAndVariantRunMarkerRunner import \
    OptimizeLFNreadCountAndVariantRunMarkerRunner
from vtam.utils.constants import get_params_default_dic


class TestOptimizeLFNreadCountAndLFNvariantRunner(unittest.TestCase):

    """Will test main commands based on a complete test dataset"""

    def setUp(self):

        nijk_path = os.path.join(os.path.dirname(__file__), 'nijk.tsv')
        known_occurrences_path = os.path.join(os.path.dirname(__file__), 'known_occurrences.tsv')

        nijk_df = pandas.read_csv(nijk_path, header=0, sep="\t")
        known_occurrences_df = pandas.read_csv(known_occurrences_path, header=0, sep="\t")

        self.nijk_df = nijk_df.loc[(nijk_df.run_id == 1) & (nijk_df.marker_id == 1)]  # one marker
        self.known_occurrences_df = known_occurrences_df.loc[
            (known_occurrences_df.run_id == 1) & (known_occurrences_df.marker_id == 1), ]  # one marker

        params_default_dic = get_params_default_dic()
        self.optimize_params_dic = {'lfn_ni_cutoff': params_default_dic['lfn_variant_cutoff'],
                                    'lfn_nik_cutoff': params_default_dic['lfn_variant_replicate_cutoff'],
                                    'lfn_njk_cutoff': params_default_dic['lfn_biosample_replicate_cutoff'],
                                    'lfn_nijk_cutoff': params_default_dic['lfn_read_count_cutoff'],
                                    'min_replicate_number': params_default_dic['min_replicate_number']}

        self.lfn_ni_cutoff = params_default_dic['lfn_variant_cutoff']
        self.lfn_nik_cutoff = params_default_dic['lfn_variant_replicate_cutoff']
        self.lfn_njk_cutoff = params_default_dic['lfn_biosample_replicate_cutoff']
        self.lfn_nijk_cutoff = params_default_dic['lfn_read_count_cutoff']
        self.min_replicate_number = params_default_dic['min_replicate_number']

        lfn_nijk_cutoff_lst = OptimizeLFNreadCountAndVariantRunMarkerRunner.get_lfn_nijk_cutoff_lst(150, 200, 3)
        lfn_ni_nik_cutoff_lst = OptimizeLFNreadCountAndVariantRunMarkerRunner.get_lfn_ni_nik_cutoff_lst(0.01, 0.51, 3)

        self.optim_run_marker_obj = OptimizeLFNreadCountAndVariantRunMarkerRunner(
            nijk_df=self.nijk_df, known_occurrences_df=self.known_occurrences_df,
            lfn_nijk_cutoff_lst=lfn_nijk_cutoff_lst, lfn_ni_nik_cutoff_lst=lfn_ni_nik_cutoff_lst)

    def test_count_keep_max_count_delete_max(self):

        count_keep_max = self.optim_run_marker_obj.get_count_keep_max()
        count_delete_max = self.optim_run_marker_obj.get_count_delete_max()

        self.assertEqual(count_keep_max, 6)
        self.assertEqual(count_delete_max, 0)

    def test_get_lst_lfn_nijk_cutoff(self):

        lfn_nijk_cutoff_lst = self.optim_run_marker_obj.get_lst_one_par_lfn_nijk_cutoff(
            lfn_ni_cutoff=self.lfn_ni_cutoff, lfn_nik_cutoff=self.lfn_nik_cutoff,
            lfn_njk_cutoff=self.lfn_njk_cutoff, lfn_nijk_cutoff=self.lfn_nijk_cutoff, min_replicate_number=self.min_replicate_number)

        self.assertEqual(lfn_nijk_cutoff_lst, [150, 170])

    def test_get_lst_one_par_lfn_ni_cutoff(self):

        lfn_ni_cutoff = self.lfn_ni_cutoff
        lfn_nik_cutoff = None

        lfn_ni_nik_cutoff_lst = self.optim_run_marker_obj.get_lst_one_par_lfn_ni_nik_cutoff(
            lfn_ni_cutoff=lfn_ni_cutoff, lfn_nik_cutoff=lfn_nik_cutoff,
            lfn_njk_cutoff=self.lfn_njk_cutoff, lfn_nijk_cutoff=self.lfn_nijk_cutoff, min_replicate_number=self.min_replicate_number)

        self.assertEqual(lfn_ni_nik_cutoff_lst, [0.01, 0.177])

    def test_get_lst_one_par_lfn_ni_nik_cutoff(self):

        lfn_ni_cutoff = None
        lfn_nik_cutoff = self.lfn_nik_cutoff

        lfn_ni_nik_cutoff_lst = self.optim_run_marker_obj.get_lst_one_par_lfn_ni_nik_cutoff(
            lfn_ni_cutoff=lfn_ni_cutoff, lfn_nik_cutoff=lfn_nik_cutoff,
            lfn_njk_cutoff=self.lfn_njk_cutoff, lfn_nijk_cutoff=self.lfn_nijk_cutoff, min_replicate_number=self.min_replicate_number)

        self.assertEqual(lfn_ni_nik_cutoff_lst, [0.01, 0.177, 0.344])

    def test_get_df_optim_lfn_readcount_variant_cutoff(self):

        lfn_ni_cutoff = self.lfn_ni_cutoff
        lfn_nik_cutoff = None

        out_two_pars_df = self.optim_run_marker_obj.get_df_optim_lfn_readcount_variant_replicate_cutoff(
            lfn_ni_cutoff=lfn_ni_cutoff, lfn_nik_cutoff=lfn_nik_cutoff,
            lfn_njk_cutoff=self.lfn_njk_cutoff, lfn_nijk_cutoff=self.lfn_nijk_cutoff, min_replicate_number=self.min_replicate_number)

        out_two_pars_df_bak = """   occurrence_nb_keep  occurrence_nb_delete  lfn_nijk_cutoff  lfn_ni_nik_cutoff
3                   6                     0              170              0.177
2                   6                     0              170              0.010
1                   6                     0              150              0.177
0                   6                     0              150              0.010"""
        self.assertEqual(out_two_pars_df.to_string(), out_two_pars_df_bak)

    def test_get_df_optim_lfn_readcount_variant_replicate_cutoff(self):

        lfn_ni_cutoff = None
        lfn_nik_cutoff = self.lfn_nik_cutoff

        out_two_pars_df = self.optim_run_marker_obj.get_df_optim_lfn_readcount_variant_replicate_cutoff(
            lfn_ni_cutoff=lfn_ni_cutoff, lfn_nik_cutoff=lfn_nik_cutoff,
            lfn_njk_cutoff=self.lfn_njk_cutoff, lfn_nijk_cutoff=self.lfn_nijk_cutoff, min_replicate_number=self.min_replicate_number)

        out_two_pars_df_bak = """   occurrence_nb_keep  occurrence_nb_delete  lfn_nijk_cutoff  lfn_ni_nik_cutoff
5                   6                     0              170              0.344
4                   6                     0              170              0.177
3                   6                     0              170              0.010
2                   6                     0              150              0.344
1                   6                     0              150              0.177
0                   6                     0              150              0.010"""
        self.assertEqual(out_two_pars_df.to_string(), out_two_pars_df_bak)
