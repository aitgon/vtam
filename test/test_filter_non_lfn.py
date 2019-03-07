import os
import pandas
from unittest import TestCase

from wopmetabarcoding.utils.PathFinder import PathFinder
from wopmetabarcoding.wrapper.FilterLFNutilities import FilterLFNRunner
from wopmetabarcoding.wrapper.FilterMinReplicateNumber import f9_delete_min_replicate_number
from wopmetabarcoding.wrapper.FilterNonLFN import FilterNonLFNRunner


class TestFilterNonLFN(TestCase):

    def setUp(self):
        pass

    def test_01_f9_delete_min_repln(self):
        #
        variant_read_count_df = pandas.DataFrame({
            'run_id': [1]*3,
            'marker_id': [1]*3,
            'variant_id': [1]*2 + [2]*1,
            'biosample_id':[1,1,1],
            'replicate_id':[1,2,1],
            'read_count':[
                10,5,249,
                  ],
        })
        min_replicate_number = 2
        #
        df_filter_output = f9_delete_min_replicate_number(variant_read_count_df = variant_read_count_df, min_replicate_number=min_replicate_number)
        #
        self.assertTrue(not df_filter_output.loc[
                            (df_filter_output.run_id == 1)
                            & (df_filter_output.marker_id == 1)
                            & (df_filter_output.variant_id == 1)
                         & (df_filter_output.biosample_id == 1)
                         & (df_filter_output.replicate_id == 1)
                         & (df_filter_output.filter_id == 9),
                                   'filter_delete'].values[0])
        self.assertTrue(df_filter_output.loc[
                            (df_filter_output.run_id == 1)
                            & (df_filter_output.marker_id == 1)
                            & (df_filter_output.variant_id == 2)
                         & (df_filter_output.biosample_id == 1)
                         & (df_filter_output.replicate_id == 1)
                         & (df_filter_output.filter_id == 9),
                                  'filter_delete'].values[0])
