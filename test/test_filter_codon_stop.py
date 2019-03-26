
from unittest import TestCase

import Bio
import pandas
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from wopmetabarcoding.wrapper.FilterCodonStop import f14_filter_chimera


class TestFilterChimera(TestCase):

    def setUp(self):
        # Input from min_replicate_number
        # Variants 1 and 2 are ok but 3-5 are chimeras
        self.variant_df = pandas.DataFrame({
            'id': list(range(1,6)),
            'sequence': [ 'AATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATATTTAGTTGG',
                          'ATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATATTTAGTTGG',
                          'AATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATAATTAGTTGG',
                          'ATAAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATAATTAGTTGG',
                          'TAAGTATTTTGACTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATAATTAGTTGG',
                         ],
        })
        #
        self.variant_read_count_df = pandas.DataFrame({
            'run_id': [1] * 12,
            'marker_id': [1] * 12,
            'variant_id': [7] * 1 + [6] * 1 + [1] * 2 + [2] * 2 + [3] * 2 + [4] * 2 + [5] * 2,
            'biosample_id': [1] * 12,
            'replicate_id': [1, 2] * 6,
            'read_count': [
                25, 25, 350, 360, 335, 325, 350, 350, 325, 325, 35, 25
            ],
        })

    def test_02_f14_codon_stop(self):
        """""
         Function searching stop codon the different reading frame of a sequence and tag/delete variant if there is
         stop codon in the 3 reading frames
         :param df_codon_stop_per_genetic_code: data frame which contains all codon stop classified by genetic code
         :param genetic_code: genetic code to search stop codon in the df_codon_stop_per_genetic_code dataframe
         :param delete_var: option which define if the variants must be deleted or tagged
         :return: void

        """
        # Input
        # self.variant_df
        # self.variant_read_count_df
        # transl_table=1
        #

        #
        df= f14_filter_chimera(self.variant_read_count_df,self.variant_df, 10)



