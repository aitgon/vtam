import os
import pathlib
import shlex
import shutil
import subprocess
import sys
import tarfile
import unittest
import urllib.request

from vtam.utils.PathManager import PathManager
from vtam.utils.constants import sorted_tar_gz_url1, sorted_tar_gz_url2, sorted_tar_gz_url3
from vtam.utils.MyProgressBar import MyProgressBar


class TestCmdFilterLFN(unittest.TestCase):

    """Will test main commands based on a complete test dataset"""

    @classmethod
    def setUpClass(cls):

        cmd = '{} -m pip install . -q --upgrade'.format(sys.executable)
        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=PathManager.get_project_path())

        cls.package_path = os.path.join(PathManager.get_package_path())
        cls.test_path = os.path.join(PathManager.get_test_path())
        cls.outdir_path = os.path.join(cls.test_path, 'outdir')
        cls.outdir_data_path = os.path.join(cls.outdir_path, 'data')
        shutil.rmtree(cls.outdir_path, ignore_errors=True)
        pathlib.Path(cls.outdir_data_path).mkdir(parents=True, exist_ok=True)

        ############################################################################################
        #
        # Download sorted reads dataset
        #
        ############################################################################################

        sorted_tar_path = os.path.join(cls.package_path, "..", "data", "sorted.tar.gz")
        # Test first in local dir, otherwise in the remote URLs
        if not os.path.isfile(sorted_tar_path) or pathlib.Path(sorted_tar_path).stat().st_size < 1000000:
            try:
                urllib.request.urlretrieve(sorted_tar_gz_url1, sorted_tar_path, MyProgressBar())
            except Exception:
                try:
                    urllib.request.urlretrieve(sorted_tar_gz_url2, sorted_tar_path, MyProgressBar())
                except Exception:
                    urllib.request.urlretrieve(sorted_tar_gz_url3, sorted_tar_path, MyProgressBar())
        tar = tarfile.open(sorted_tar_path, "r:gz")
        tar.extractall(path=cls.outdir_data_path)
        tar.close()

    def setUp(self):
        self.outdir_thistest_path = os.path.join(self.outdir_path, 'thistest')
        # during development of the test, this prevents errors
        pathlib.Path(self.outdir_thistest_path).mkdir(parents=True, exist_ok=True)
        os.environ['VTAM_LOG_VERBOSITY'] = str(10)

        ############################################################################################
        #
        # Paths
        #
        ############################################################################################

        self.asvtable_path = os.path.join(self.outdir_path, "asvtable_default.tsv")

        self.args = {}
        self.args['sortedinfo'] = os.path.join(os.path.dirname(__file__), "sortedinfo.tsv")
        self.args['sorteddir'] = os.path.join(self.outdir_data_path, 'sorted')
        self.args['optimize_lfn_variant_specific'] = os.path.join(
            self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/optimize_lfn_variant_specific.tsv")
        self.args['optimize_lfn_variant_replicate_specific'] = os.path.join(
            self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/optimize_lfn_variant_replicate_specific.tsv")
        self.args['params_lfn_variant'] = os.path.join(os.path.dirname(__file__), "params_lfn_variant.yml")
        self.args['params_lfn_variant_replicate'] = os.path.join(os.path.dirname(__file__), "params_lfn_variant_replicate.yml")

    def test_filter_params_lfn_variant_replicate_fails1(self):

        ############################################################################################
        #
        # fails
        #
        ############################################################################################

        cmd = "vtam filter --db db.sqlite --sortedinfo {sortedinfo} --sorteddir {sorteddir} " \
              "--asvtable asvtable_default.tsv --params {params_lfn_variant} --until FilterLFN " \
              "--lfn_variant_replicate".format(**self.args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        result = subprocess.run(args=args, cwd=self.outdir_thistest_path)

        self.assertEqual(result.returncode, 1)

    def test_filter_params_lfn_variant_replicate_fails2(self):

        ############################################################################################
        #
        # fails
        #
        ############################################################################################

        cmd = "vtam filter --db db.sqlite --sortedinfo {sortedinfo} --sorteddir {sorteddir} " \
              "--asvtable asvtable_default.tsv --params {params_lfn_variant_replicate} --until FilterLFN " \
              "".format(**self.args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        result = subprocess.run(args=args, cwd=self.outdir_thistest_path)

        self.assertEqual(result.returncode, 1)

    def test_filter_params_lfn_variant_replicate_succeeds1(self):

        ############################################################################################
        #
        # succeeds
        #
        ############################################################################################

        cmd = "vtam filter --db db.sqlite --sortedinfo {sortedinfo} --sorteddir {sorteddir} " \
              "--asvtable asvtable_default.tsv --params {params_lfn_variant} --until FilterLFN " \
              "".format(**self.args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        result = subprocess.run(args=args, cwd=self.outdir_thistest_path)

        self.assertEqual(result.returncode, 0)

    def test_filter_params_lfn_variant_replicate_succeeds2(self):

        ############################################################################################
        #
        # succeeds
        #
        ############################################################################################

        cmd = "vtam filter --db db.sqlite --sortedinfo {sortedinfo} --sorteddir {sorteddir} " \
              "--asvtable asvtable_default.tsv --params {params_lfn_variant_replicate} --until FilterLFN " \
              "--lfn_variant_replicate".format(**self.args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        result = subprocess.run(args=args, cwd=self.outdir_thistest_path)

        self.assertEqual(result.returncode, 0)

    def test_filter_fasta_over_two_lines(self):

        cmd = "vtam filter --db db.sqlite --sortedinfo {sortedinfo} --sorteddir {sorteddir} " \
              "--asvtable asvtable_default.tsv --params {params_lfn_variant_replicate} --until FilterLFN " \
              "--lfn_variant_replicate".format(**self.args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        result = subprocess.run(args=args, cwd=self.outdir_thistest_path)
        self.assertEqual(result.returncode, 0)

    def tearDown(self):

        shutil.rmtree(self.outdir_thistest_path, ignore_errors=True)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.outdir_path, ignore_errors=True)
