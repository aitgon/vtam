import multiprocessing
import os
import pathlib

import pandas

from unittest import TestCase
from vtam.utils.PathManager import PathManager
from vtam.utils.FilterPCRerrorRunner import FilterPCRerrorRunner

class TestFilterPCRError(TestCase):

    def setUp(self):
        os.environ['VTAM_THREADS'] = str(multiprocessing.cpu_count())

        # Input from min_replicate_number
        self.variant_expected_df = pandas.DataFrame({
            'sequence': [
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTG',
                         ],
        }, index=[1])
        # Input from min_replicate_number
        self.variant_unexpected_df = pandas.DataFrame({
            'sequence': [
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTT',
                         ],
        }, index=[2])
        #
        self.variant_read_count_df = pandas.DataFrame({
            'run_id': [1]*4,
            'marker_id': [1]*4,
            'biosample_id': [1]*4,
            'replicate_id': [1, 2]*2,
            'variant_id': [1]*2 + [2]*2,
            'read_count':[
                100, 100, 1, 1,
                  ],
        })

        self.this_tempdir = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__))
        pathlib.Path(self.this_tempdir).mkdir(parents=True, exist_ok=True)

    def test_get_vsearch_alignement_df(self):
        filter_pcr_error_runner = FilterPCRerrorRunner(variant_expected_df=self.variant_expected_df, variant_unexpected_df=self.variant_unexpected_df,
                                                       variant_read_count_df=self.variant_read_count_df, tmp_dir=self.this_tempdir)
        vsearch_alignement_df = filter_pcr_error_runner.get_vsearch_alignement_df()

        vsearch_alignement_df_bak_str = """   variant_id_unexpected  variant_id_expected  alnlen  ids  mism  gaps
0                      2                    1      32   31     1     0"""
        self.assertTrue(vsearch_alignement_df_bak_str==vsearch_alignement_df.to_string())

    def test_get_filter_pcr_error_df(self):

        filter_pcr_error_runner = FilterPCRerrorRunner(variant_expected_df=self.variant_expected_df, variant_unexpected_df=self.variant_unexpected_df,
                                                       variant_read_count_df=self.variant_read_count_df, tmp_dir=self.this_tempdir)
        variant_unexpected_to_expected_ratio_df = filter_pcr_error_runner.get_variant_unexpected_to_expected_ratio_df()

        variant_unexpected_to_expected_ratio_df_bak = pandas.DataFrame({
            'run_id' : [1],
            'marker_id' : [1],
            'biosample_id' : [1],
            'variant_id_expected' : [1],
            'N_ij_expected' : [200],
            'variant_id_unexpected' : [2],
            'N_ij_unexpected' : [2],
            'N_ij_unexpected_to_expected_ratio' : [0.01],
        })
        #
        pandas.testing.assert_frame_equal(variant_unexpected_to_expected_ratio_df, variant_unexpected_to_expected_ratio_df_bak)



    def test_get_filter_output_df(self):

        filter_pcr_error_runner = FilterPCRerrorRunner(variant_expected_df=self.variant_expected_df, variant_unexpected_df=self.variant_unexpected_df,
                                                       variant_read_count_df=self.variant_read_count_df, tmp_dir=self.this_tempdir)
        pcr_error_var_prop = 0.05

        filter_output_df = filter_pcr_error_runner.get_filter_output_df(pcr_error_var_prop)

        filter_output_df_bak = pandas.DataFrame({
            'run_id': [1]*4,
            'marker_id': [1]*4,
            'biosample_id': [1]*4,
            'replicate_id': [1, 2]*2,
            'variant_id': [1]*2 + [2]*2,
            'read_count': [100, 100, 1, 1,],
            'filter_delete': [False, False, True, True],
        })
        #
        pandas.testing.assert_frame_equal(filter_output_df, filter_output_df_bak)

