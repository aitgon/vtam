import os
import pandas
import unittest

from vtam.utils.PathManager import PathManager
from vtam.utils.SequenceClusterer import SequenceClusterer


class TestSequenceClusterer(unittest.TestCase):

    def setUp(self):

        self.test_path = os.path.join(PathManager.get_test_path())

    def test_compute_clusters_097(self):

        asvtable_default_tsv = os.path.join(self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/asvtable_default.tsv")

        asvtable_df = pandas.read_csv(asvtable_default_tsv, sep='\t')
        asvtable_df.rename({'variant': 'variant_id'}, axis=1, inplace=True)

        cluster_count_097_df = SequenceClusterer(asvtable_df, cluster_identity=0.97).compute_clusters()
        cluster_count_095_df = SequenceClusterer(asvtable_df, cluster_identity=0.95).compute_clusters()

        clusterid_097_list = [121, 1234, 1234, 1234, 1234, 1241, 1309, 1309, 1328, 144, 144, 1371, 4347, 4347, 1651, 1651, 2031, 2031, 1959, 2124, 2138, 2138, 2371, 2371, 2371, 2371, 2524, 2542, 2684, 2684, 3679, 3679, 2790, 2801, 2801, 2816, 3101, 3101, 3101, 32, 3299, 3307, 4420, 4420, 380, 380, 400, 5031, 5031, 4073, 4134, 422, 4265, 4265, 4355, 4378, 4445, 456, 456, 456, 474, 493, 507, 510, 521, 538, 553, 58, 603, 954]
        self.assertEqual(cluster_count_097_df.clusterid.tolist(), clusterid_097_list)
        self.assertNotEqual(cluster_count_095_df.clusterid.tolist(), clusterid_097_list)

        clustersize_097_list = [1, 4, 4, 4, 4, 1, 2, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1, 1, 2, 2, 4,
         4, 4, 4, 1, 1, 2, 2, 2, 2, 1, 2, 2, 1, 3, 3, 3, 1, 1, 1, 2, 2, 2, 2,
         1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1]
        self.assertEqual(cluster_count_097_df.clustersize.tolist(), clustersize_097_list)
        self.assertNotEqual(cluster_count_095_df.clustersize.tolist(), clustersize_097_list)
