
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


        number_codon_stop_1 = 0
        number_codon_stop_2 = 0
        number_codon_stop_3 = 0
        stop_codons = ["TAG", "TAA"]
        sequence = self.variant_df.sequence[0]

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
        import pdb;
        pdb.set_trace()




