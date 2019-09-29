
from unittest import TestCase
from wopmetabarcoding.utils.PathManager import PathManager

import os
from wopmetabarcoding.utils.utilities import tempdir
import pandas

from wopmetabarcoding.wrapper.FilterIndel import f13_filter_indel


class TestIndel(TestCase):

    def setUp(self):
        # Input from min_replicate_number
        # Variants 1 and 2 are ok but 3-5 are chimeras
        self.variant_df = pandas.DataFrame({
            'id': list(range(1,8)),
            'sequence': [
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTAGCTGTAACTTTGTCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTA',
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTAGATTACTATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTA',
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTATTTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTAA',
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTA',
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTTTCAATTTTATAACACTATTGTAACTGCTCATGGATTATTA',
                'TGTTCTTTATTTATTATTGCTGGTTTTGCTGGTGTTTTAGCTGTAACTTTATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTATCAATTTTATAACACTATTGTAACTGCTCATGGATTATTA',
                'TGTTCTTTATTTATTATTTGCTGGTTTTGCTGGTGTTTTAGCTGATCATTATTAATTAGATTACAATTAGTTGCTACTGGGTATGGATGATTAGCTTTGAATTTTCAATTTTATAACACTATTGTAACTGCTCATGGATTATTA'

                         ],
        })
        # output
        # seq_indel1(var1) and seq_indel7(var7) should be detected as potential pseudogenes/sequencing errors
        #
        self.variant_read_count_df = pandas.DataFrame({
            'run_id': [1] * 12,
            'marker_id': [1]*12,
            'variant_id': [7]*1 + [6]*1 + [1]*2 + [2]*2 + [3]*2 + [4]*2 + [5]*2,
            'biosample_id': [1] * 12,
            'replicate_id': [1, 2] * 6,
            'read_count':[
                25, 25, 350, 360, 335, 325, 350, 350, 325, 325, 35, 25
                  ],
        })

        self.tempdir = os.path.join(tempdir, "FilterUtilities", self.__class__.__name__)
        PathManager.mkdir_p(self.tempdir)


    def test_01_f13_indel(self):

        df_out = f13_filter_indel(self.variant_read_count_df, self.variant_df)
        # df_out1 = f13_filter_indel(self.variant_read_count_df, self.variant_df1)
        # Variant 1 passes
        self.assertTrue(not df_out.loc[(df_out.run_id == 1)
                                         & (df_out.marker_id == 1)
                                         & (df_out.variant_id == 1)
                                         & (df_out.biosample_id == 1)
                                         & (df_out.replicate_id == 1),
                                         # & (df_out.filter_id == 13),
                                         'filter_delete'].values[0])
        # Variant 7 passes
        self.assertTrue(df_out.loc[(df_out.run_id == 1)
                                         & (df_out.marker_id == 1)
                                         & (df_out.variant_id == 7)
                                         & (df_out.biosample_id == 1)
                                         & (df_out.replicate_id == 1),
                                         # & (df_out.filter_id == 13),
                                         'filter_delete'].values[0])








