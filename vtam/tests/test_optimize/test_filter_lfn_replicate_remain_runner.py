import os
import pandas
import unittest

from vtam.utils.constants import get_params_default_dic
from vtam.utils.FilterLFNreplicateRemainRunner import FilterLFNreplicateRemainRunner


class TestFilterLFNreplicateRemainRunner(unittest.TestCase):

    """Will test main commands based on a complete test dataset"""

    def setUp(self):

        nijk_path = os.path.join(os.path.dirname(__file__), 'nijk.tsv')
        known_occurrences_path = os.path.join(os.path.dirname(__file__), 'known_occurrences.tsv')

        nijk_df = pandas.read_csv(nijk_path, header=0, sep="\t")
        known_occurrences_df = pandas.read_csv(known_occurrences_path, header=0, sep="\t")

        self.nijk_df = nijk_df.loc[(nijk_df.run_id == 1) & (nijk_df.marker_id == 1)]  # one marker
        self.known_occurrs_run_marker_df = known_occurrences_df.loc[
            (known_occurrences_df.run_id == 1) & (known_occurrences_df.marker_id == 1), ]  # one marker

        params_default_dic = get_params_default_dic()
        self.optimize_params_dic = {'lfn_ni_cutoff': params_default_dic['lfn_variant_cutoff'],
                                    'lfn_nik_cutoff': params_default_dic['lfn_variant_replicate_cutoff'],
                                    'lfn_njk_cutoff': params_default_dic['lfn_sample_replicate_cutoff'],
                                    'lfn_nijk_cutoff': params_default_dic['lfn_read_count_cutoff'],
                                    'min_replicate_number': params_default_dic['min_replicate_number']}

    def test_lfn_variant(self):

        lfn_ni_cutoff_lst = [0.2, 0.3, 0.35, 0.4]
        self.optimize_params_dic['lfn_nik_cutoff'] = None
        self.optimize_params_dic['lfn_nijk_cutoff'] = 100
        count_keep_lst = [6, 5, 2, 0]

        for i, lfn_ni_cutoff in enumerate(lfn_ni_cutoff_lst):
            self.optimize_params_dic['lfn_ni_cutoff'] = lfn_ni_cutoff
            count_keep, count_delete = FilterLFNreplicateRemainRunner(nijk_df=self.nijk_df, **self.optimize_params_dic)\
                .count_keep_delete(known_occurrences_df=self.known_occurrs_run_marker_df)
            self.assertEqual(count_keep, count_keep_lst[i])

    def test_lfn_variant_replicate(self):

        lfn_nik_cutoff_lst = [x*3 for x in [0.2, 0.3, 0.35, 0.4]]
        self.optimize_params_dic['lfn_ni_cutoff'] = None
        self.optimize_params_dic['lfn_nijk_cutoff'] = 100
        count_keep_lst = [6, 6, 0, 0]

        for i, lfn_nik_cutoff in enumerate(lfn_nik_cutoff_lst):
            self.optimize_params_dic['lfn_nik_cutoff'] = lfn_nik_cutoff
            count_keep, count_delete = FilterLFNreplicateRemainRunner(nijk_df=self.nijk_df, **self.optimize_params_dic)\
                .count_keep_delete(known_occurrences_df=self.known_occurrs_run_marker_df)
            self.assertEqual(count_keep, count_keep_lst[i])
