import os
import pandas
from unittest import TestCase

from vtam.utils.FilterPCRerrorRunner import FilterPCRerrorRunner
from vtam.utils.PathManager import PathManager
# from vtam.wrapper.FilterPCRerror import f10_pcr_error_analyze_vsearch_output_df, f10_pcr_error_run_vsearch


class TestFilterPCRError(TestCase):

    def setUp(self):
        # Input from min_replicate_number
        self.variant_df = pandas.DataFrame({
            'id': list(range(1,5)),
            'sequence': [
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTAATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATTATTAGTTGG',
                'TGTTCTTTATTTATTATTTGATGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTAATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATTATTAGTTGG',
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTCGCTGTAACTTTATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTAATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATTATTAGTTGG',
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTCGCTGTAACTTTATCATTATCAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTAATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATTATTAGTTGG',
                         ],
        })
        #
        self.variant_read_count_df = pandas.DataFrame({
            'run_id': [1]*8,
            'marker_id': [1]*8,
            'biosample_id': [1]*8,
            'replicate_id': [1, 2]*4,
            'variant_id': [1]*2 + [2]*2 + [3]*2 + [4]*2,
            'read_count':[
                350, 300, 300, 220, 60, 0, 2, 0,
                  ],
        })
        #
        # Output of vsearch
        self.pcr_error_vsearch_output_df = pandas.DataFrame({
            'query' : [3, 3, 3, 3, 2, 2, 2, 2, 4, 4, 4, 4, 1, 1],
            'target' : [3, 1, 4, 2, 2, 1, 3, 4, 4, 3, 1, 2, 1, 2],
            'alnlen' : [300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300],
            'ids' : [300, 299, 299, 298, 300, 299, 298, 297, 300, 299, 298, 297, 300, 299],
            'mism' : [0, 1, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1],
            'gaps' : [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        })

        self.this_step_tmp_dir = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__))
        PathManager.mkdir_p(self.this_step_tmp_dir)



    def test_get_vsearch_alignement_df(self):

        filter_pcr_error_runner = FilterPCRerrorRunner(variant_expected_df=self.variant_df, variant_unexpected_df=self.variant_df,
                             variant_read_count_df=self.variant_read_count_df, tmp_dir=self.this_step_tmp_dir)
        vsearch_alignement_df = filter_pcr_error_runner.get_vsearch_alignement_df()
        self.assertTrue(sorted(vsearch_alignement_df.ids.unique().tolist()) == [297, 298, 299, 300])

    def test_get_filter_output_df(self):

        filter_pcr_error_runner = FilterPCRerrorRunner(variant_expected_df=self.variant_df, variant_unexpected_df=self.variant_df,
                             variant_read_count_df=self.variant_read_count_df, tmp_dir=self.this_step_tmp_dir)
        #
        pcr_error_var_prop = 0.05
        filter_output_df = filter_pcr_error_runner.get_filter_output_df(pcr_error_var_prop)
        self.assertTrue(filter_output_df.read_count.tolist() == [350, 300, 300, 220, 60, 0, 2, 0])
        self.assertTrue(filter_output_df.filter_delete.tolist() == [False, False, False, False, False, False, True, True])
