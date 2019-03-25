
from unittest import TestCase

import pandas

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

        number_codon_stop_1 = 0
        number_codon_stop_2 = 0
        number_codon_stop_3 = 0
        df = self.variant_df.copy()
        stop_codons = ["TAG","TAA"]
        df['number_codon_stop_1'] = 0
        df['number_codon_stop_2'] = 0
        df['number_codon_stop_3'] = 0
        df['orf_ref']= 'none'
        for row in df.iterrows():
            sequence = row[1].sequence
            id = row[1].id

            for i in range(0, len(sequence), 3):
                codon = sequence[i:i + 3]
                if codon in stop_codons:
                  number_codon_stop_1 = number_codon_stop_1 + 1
            for i in range(1, len(sequence), 3):
                codon = sequence[i:i + 3]
                if codon in stop_codons:
                   number_codon_stop_2 = number_codon_stop_2 + 1
            for i in range(2, len(sequence), 3):
                codon = sequence[i:i + 3]
                if codon in stop_codons:
                  number_codon_stop_3 = number_codon_stop_3 + 1
            df.loc[df.id==id,'number_codon_stop_1']=number_codon_stop_1
            df.loc[df.id == id, 'number_codon_stop_2'] = number_codon_stop_2
            df.loc[df.id == id, 'number_codon_stop_3'] = number_codon_stop_3
            #
            m = min(number_codon_stop_1,number_codon_stop_2)

            if (m==number_codon_stop_1) &( m==0):
                df.loc[df.id==id,'orf_ref']='orf_1'
            if (m==number_codon_stop_2) & (m==0):
                df.loc[df.id==id,'orf_ref']='orf_2'
            if (m==number_codon_stop_3) & (m==0):
                df.loc[df.id == id, 'orf_ref'] = 'orf_3'
            #
            number_codon_stop_1 = 0
            number_codon_stop_2 = 0
            number_codon_stop_3 = 0
        import pdb;
        pdb.set_trace()


