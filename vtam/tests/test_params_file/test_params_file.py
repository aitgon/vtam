import os
import unittest

from vtam.utils.ParamsFile import ParamsFile
from vtam.utils.PathManager import PathManager


class TestParamsFile(unittest.TestCase):

    def setUp(self):

        package_path = PathManager.get_package_path()
        self.params = os.path.join(package_path, "doc/data/params_mfzr.yml")
        self.params_wrong = os.path.join(os.path.dirname(__file__), "params_wrong.yml")

    def test_is_valid(self):

        self.assertTrue(ParamsFile(self.params).is_valid())

    def test_params_is_not_valid(self):

        with self.assertRaises(SystemExit):
            ParamsFile(self.params_wrong).is_valid()

    def test_get_params_dic(self):

        self.assertTrue(ParamsFile(self.params).get_params_dic()['lfn_biosample_replicate_cutoff']
                        == 0.003)
