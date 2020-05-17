import pandas
import unittest

from vtam.utils.FilterIndelRunner import FilterIndelRunner


class TestIndel(unittest.TestCase):

    def setUp(self):

        self.variant_df = pandas.DataFrame(
            data={
                'sequence': [
                    'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTAGCTGTAACTTTGTCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTA',
                    'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTAGATTACTATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTA',
                    'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTATTTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTAA',
                    'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTA',
                    'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTTTCAATTTTATAACACTATTGTAACTGCTCATGGATTATTA',
                    'TGTTCTTTATTTATTATTGCTGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTA',
                    'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTAGCTGATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTTTCAATTTTATAACACTATTGTAACTGCTCATGGATTATTA'],
            },
            index=list(
                range(
                    1,
                    8)))
        #
        self.variant_read_count_df = pandas.DataFrame({
            'run_id': [1] * 7,
            'marker_id': [1] * 7,
            'variant_id': list(range(1, 8)),
            'biosample_id': [1] * 7,
            'replicate': [1] * 7,
            'read_count': [
                25, 25, 350, 360, 335, 325, 350
            ],
        })

    def test_01(self):

        variant_read_count_delete_df = FilterIndelRunner(
            self.variant_read_count_df).get_variant_read_count_delete_df(
            variant_df=self.variant_df, skip_filter_indel=0)

        # Variant 1 passes
        self.assertFalse(
            variant_read_count_delete_df.loc[
                (variant_read_count_delete_df.run_id == 1) & (
                    variant_read_count_delete_df.marker_id == 1) & (
                    variant_read_count_delete_df.variant_id == 1) & (
                    variant_read_count_delete_df.biosample_id == 1) & (
                        variant_read_count_delete_df.replicate == 1),
                'filter_delete'].values[0])
        # Variant 7 passes
        self.assertTrue(
            variant_read_count_delete_df.loc[
                (variant_read_count_delete_df.run_id == 1) & (
                    variant_read_count_delete_df.marker_id == 1) & (
                    variant_read_count_delete_df.variant_id == 7) & (
                    variant_read_count_delete_df.biosample_id == 1) & (
                        variant_read_count_delete_df.replicate == 1),
                'filter_delete'].values[0])
