import pathlib
from unittest import TestCase

from vtam.utils.PathManager import PathManager
import os
import pandas
from vtam.wrapper.FilterRenkonen import f12_filter_delete_renkonen

from vtam.wrapper.FilterRenkonen import renkonen_distance


class TestFilterRenkonen(TestCase):

    def setUp(self):
        self.variant_read_count_df = pandas.DataFrame({
            'run_id': [1] * 6,
            'marker_id': [1] * 6,
            'variant_id': [1] * 3 + [2] * 3,
            'biosample_id': [1] * 6,
            'replicate': [1, 2, 3] * 2,
            'read_count': [
                5180, 5254, 9378, 193, 99, 209
            ],
        })

        self.this_tempdir = os.path.join(PathManager.instance().get_tempdir(), self.__class__.__name__)
        pathlib.Path(self.this_tempdir).mkdir(exist_ok=True)


    def test_f12_delete_filter_renkonen(self):
        #
        Rthr = 0.005
        filter_output_df = f12_filter_delete_renkonen(self.variant_read_count_df, Rthr)
        self.assertTrue(filter_output_df.loc[(filter_output_df.run_id == 1)
                                             & (filter_output_df.marker_id == 1)
                                             & (filter_output_df.variant_id == 1)
                                             & (filter_output_df.biosample_id == 1)
                                             & (filter_output_df.replicate == 1),
                                             # & (filter_output_df.filter_id == 12),
                                             'filter_delete'].values[0])
        #
        self.assertTrue(not filter_output_df.loc[(filter_output_df.run_id == 1)
                                             & (filter_output_df.marker_id == 1)
                                             & (filter_output_df.variant_id == 2)
                                             & (filter_output_df.biosample_id == 1)
                                             & (filter_output_df.replicate == 3),
                                             # & (filter_output_df.filter_id == 12),
                                             'filter_delete'].values[0])

        # import pdb;
        # pdb.set_trace()

    def test_renkonen_distance(self):
        # Output
        # biosample_id 1, replicate 1, replicate 2, renkonen_similarity and distance 0.982573959807409 and 0.017426040192591
        # biosample_id 1, replicate 1, replicate 3, renkonen_similarity and distance 0.985880012193912 and 0.014119987806088
        run_id = 1
        marker_id = 1
        biosample_id = 1
        left_replicate = 1
        right_replicate = 2
        #
        distance_left_right = renkonen_distance(self.variant_read_count_df,run_id,marker_id,biosample_id,left_replicate,right_replicate)
        self.assertAlmostEqual(distance_left_right, 0.017426040192591)
        right_replicate = 3
        distance_left_right= renkonen_distance(self.variant_read_count_df,run_id,marker_id,biosample_id,left_replicate,right_replicate)
        self.assertAlmostEqual(distance_left_right, 0.014119987806088)

