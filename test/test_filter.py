import filecmp
import os
import pickle
import shutil
import sqlite3
from unittest import TestCase

import pandas

from wopmetabarcoding.wrapper.FilterUtilities import FilterRunner


class TestFilter(TestCase):
    def setUp(self):
        self.variant_df = pandas.DataFrame({
            'id':[1,22],
            'sequence_':["tata", "tgtg"],
        })
        self.variant_read_count_df = pandas.DataFrame({
            'variant_id':[1,1,1,1,1,1,22,22,22,22,22,22],
            'biosample_id':[1,1,1,2,2,2,1,1,1,2,2,2],
            'replicate_id':[1,2,3,1,2,3,1,2,3,1,2,3],
            'read_count':[10,5,0,249,58,185,25,58,23,10980,8999,13814],
        })
        self.marker_id = 1
        #
        self.filter_runner = FilterRunner(self.variant_df, self.variant_read_count_df, self.marker_id)

    def test_02_f2_lfn2_per_variant_mekdad(self):
        lfn_var_threshold = 0.001
        self.filter_runner.f2_lfn2_per_variant_mekdad(lfn_var_threshold)
        self.assertTrue(self.filter_runner.passed_variant_mekdad_df.f2_lfn2_per_variant_mekdad.values.tolist() == [True, True, False, True, True, True, False, True, False, True, True, True])
        self.assertFalse(self.filter_runner.passed_variant_mekdad_df.f2_lfn2_per_variant_mekdad.values.tolist() == [True, True, False, True, True, True, False, True, False, True, True, False])
