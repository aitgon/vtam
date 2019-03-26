
from unittest import TestCase

import Bio
import pandas
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq


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
        orf1_codon_stop_nb = 0
        orf2_codon_stop_nb = 0
        orf3_codon_stop_nb = 0
        #
        df = self.variant_df.copy()
        #
        variant_sequence = "AATAGTATTTTTTCTCCTTATGCCTGCTTTAATAGGTGGTTTTGGTAATTGAATAGTTCCTGTTCTAATTGGTTCTATTGATATGGCTTACCCTAGATTAAATAATATTAGTTTTTGATTATTGCCCCCTAGTTTATTATATTTAGTTGG"
        orf_frame_index = 1 # 1-based
        str(Seq(variant_sequence[orf_frame_index-1:], IUPAC.unambiguous_dna).translate(Bio.Data.CodonTable.generic_by_id[1])).count('*')
        #
        # stop_codons = ["TAG","TAA"] # Read from data/genetic_codes.csv
        # #
        # df['orf1_codon_stop_nb'] = 0
        # df['orf2_codon_stop_nb'] = 0
        # df['orf3_codon_stop_nb'] = 0
        # #
        # df['orf_index_with_0_codon_stop']= 'none'
        # #
        # import pdb; pdb.set_trace()
        # for row in df.iterrows():
        #     sequence = row[1].sequence
        #     id = row[1].id
        #
        #     for i in range(0, len(sequence), 3):
        #         codon = sequence[i:i + 3]
        #         if codon in stop_codons:
        #           orf1_codon_stop_nb = orf1_codon_stop_nb + 1
        #     for i in range(1, len(sequence), 3):
        #         codon = sequence[i:i + 3]
        #         if codon in stop_codons:
        #            orf2_codon_stop_nb = orf2_codon_stop_nb + 1
        #     for i in range(2, len(sequence), 3):
        #         codon = sequence[i:i + 3]
        #         if codon in stop_codons:
        #           orf3_codon_stop_nb = orf3_codon_stop_nb + 1
        #     df.loc[df.id==id,'orf1_codon_stop_nb']=orf1_codon_stop_nb
        #     df.loc[df.id == id, 'orf2_codon_stop_nb'] = orf2_codon_stop_nb
        #     df.loc[df.id == id, 'orf3_codon_stop_nb'] = orf3_codon_stop_nb
        #     #
        #     m = min(orf1_codon_stop_nb,orf2_codon_stop_nb)
        #
        #     if (m==orf1_codon_stop_nb) &( m==0):
        #         df.loc[df.id==id,'orf_index_with_0_codon_stop']='orf_1'
        #     if (m==orf2_codon_stop_nb) & (m==0):
        #         df.loc[df.id==id,'orf_index_with_0_codon_stop']='orf_2'
        #     if (m==orf3_codon_stop_nb) & (m==0):
        #         df.loc[df.id == id, 'orf_index_with_0_codon_stop'] = 'orf_3'
        #     #
        #     orf1_codon_stop_nb = 0
        #     orf2_codon_stop_nb = 0
        #     orf3_codon_stop_nb = 0
        # import pdb;
        # pdb.set_trace()


