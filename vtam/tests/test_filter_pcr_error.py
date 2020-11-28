import multiprocessing
import os
import pathlib

import pandas
from unittest import TestCase

from vtam.utils.RunnerFilterPCRerror import RunnerFilterPCRerror
from vtam.utils.PathManager import PathManager


class TestFilterPCRError(TestCase):

    def setUp(self):
        os.environ['VTAM_THREADS'] = str(multiprocessing.cpu_count())

        # Input from min_replicate_number
        self.variant_df = pandas.DataFrame({
            'sequence': [
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTAATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATTATTAGTTGG',
                'TGTTCTTTATTTATTATTTGATGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTAATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATTATTAGTTGG',
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTCGCTGTAACTTTATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTAATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATTATTAGTTGG',
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTCGCTGTAACTTTATCATTATCAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTAATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATTATTAGTTGG',
            ],
        }, index=list(range(1, 5)))
        #
        self.variant_read_count_df = pandas.DataFrame({
            'run_id': [1] * 8,
            'marker_id': [1] * 8,
            'sample_id': [1] * 8,
            'replicate': [1, 2] * 4,
            'variant_id': [1] * 2 + [2] * 2 + [3] * 2 + [4] * 2,
            'read_count': [
                350, 300, 300, 220, 60, 0, 2, 0,
            ],
        })

        self.this_tempdir = os.path.join(
            PathManager.instance().get_tempdir(),
            os.path.basename(__file__))
        pathlib.Path(self.this_tempdir).mkdir(parents=True, exist_ok=True)

    def test_get_vsearch_alignement_df(self):

        filter_pcr_error_runner = RunnerFilterPCRerror(
            variant_expected_df=self.variant_df,
            variant_unexpected_df=self.variant_df,
            variant_read_count_df=self.variant_read_count_df)
        vsearch_alignement_df = filter_pcr_error_runner.get_vsearch_alignement_df()
        self.assertTrue(
            sorted(
                vsearch_alignement_df.ids.unique().tolist()) == [
                297,
                298,
                299,
                300])

    def test_get_filter_output_df(self):

        filter_pcr_error_runner = RunnerFilterPCRerror(
            variant_expected_df=self.variant_df,
            variant_unexpected_df=self.variant_df,
            variant_read_count_df=self.variant_read_count_df)
        #
        pcr_error_var_prop = 0.05
        filter_output_df = filter_pcr_error_runner.get_variant_read_count_delete_df(
            pcr_error_var_prop)
        filter_output_bak_dic = {'filter_delete': {0: False,
                   1: False,
                   2: False,
                   3: False,
                   4: False,
                   5: False,
                   6: True,
                   7: True},
 'marker_id': {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1},
 'read_count': {0: 350, 1: 300, 2: 300, 3: 220, 4: 60, 5: 0, 6: 2, 7: 0},
 'replicate': {0: 1, 1: 2, 2: 1, 3: 2, 4: 1, 5: 2, 6: 1, 7: 2},
 'run_id': {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1},
 'sample_id': {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1},
 'variant_id': {0: 1, 1: 1, 2: 2, 3: 2, 4: 3, 5: 3, 6: 4, 7: 4}}
        self.assertTrue(filter_output_bak_dic ==
                        filter_output_df.to_dict())
