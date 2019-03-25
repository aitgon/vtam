
from unittest import TestCase

import pandas

from wopmetabarcoding.wrapper.FilterConsensus import f15_filter_consensus


class TestFilterConsensus(TestCase):

    def setUp(self):
        # Input from min_replicate_number
        #
        self.variant_read_count_df = pandas.DataFrame({
            'run_id': [1] * 5,
            'marker_id': [1] * 5,
            'variant_id': [1] * 2 + [2] * 3,
            'biosample_id': [1] * 2 + [2] * 2 + [1],
            'replicate_id': [1, 2] + [1, 2] + [3],
            'read_count':[
                156, 341, 99, 140, 116
                  ],
         })

    def test_02_f15_consensus(self):
        filter_output_df = f15_filter_consensus(self.variant_read_count_df)

        #
        self.assertTrue((filter_output_df.loc[(filter_output_df.run_id == 1)
                                                                         & (filter_output_df.marker_id == 1)
                                                                         & (filter_output_df.variant_id == 1)
                                                                         & (filter_output_df.biosample_id == 1),
                                                                        'read_count'].values[0]) == 248.5)
        #
        self.assertTrue((filter_output_df.loc[(filter_output_df.run_id == 1)
                                              & (filter_output_df.marker_id == 1)
                                              & (filter_output_df.variant_id == 2)
                                              & (filter_output_df.biosample_id == 2),
                                              'read_count'].values[0]) == 119.5)
         #

        self.assertTrue((filter_output_df.loc[(filter_output_df.run_id == 1)
                                              & (filter_output_df.marker_id == 1)
                                              & (filter_output_df.variant_id == 2)
                                              & (filter_output_df.biosample_id == 1),
                                              'read_count'].values[0]) == 116.0)

