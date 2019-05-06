# -*- coding: utf-8 -*-

from unittest import TestCase

class TestMakeTableOTU(TestCase):

    def setUp(self):
        pass


    def test_f01_make_table_out(self):
        # Input
        variant_dic = {
            'variant_id' : [102],
            'variant_sequence' : ['ACTTTATTTCATTTTCGGAACATTTGCAGGAGTTGTAGGAACTTTACTTTCATTATTTATTCGACTAGAATTAGCTTATCCAGGAAATCAATTTTTTTTAGGAAATCACCAACTTTATAATGTGGTTGTGACAGCACATGCTTTTATCATGATTTTTTTCATGGTTATGCCGATTTTAATC']
        }
        filter_codon_stop_dic = {'run_id': [1, 1, 1, 1, 1, 1], 'marker_id': [1, 1, 1, 1, 1, 1],
                                 'variant_id': [102, 102, 102, 102, 102, 102], 'biosample_id': [1, 1, 1, 2, 2, 2],
                                 'replicate_id': [1, 2, 3, 1, 2, 3], 'read_count': [217, 174, 150, 266, 84, 96],
                                 'fiter_id': [14, 14, 14, 14, 14, 14], 'filter_delete': [0, 0, 0, 0, 0, 0]}
        #

