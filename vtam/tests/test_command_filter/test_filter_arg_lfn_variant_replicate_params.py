import filecmp
import os
import pathlib
import shlex
import shutil
import sqlite3
import subprocess
import sys
import tarfile
import unittest
import urllib

from vtam.utils.PathManager import PathManager
from vtam.utils.constants import sorted_tar_gz_url
from urllib import request

@unittest.skipIf(request.urlopen(sorted_tar_gz_url).getcode() != 200,
                 "This test requires an internet connection!")
class TestCommands(unittest.TestCase):

    """Will test main commands based on a complete test dataset"""

    @classmethod
    def setUpClass(cls):

        # vtam needs to be in the tsv_path
        subprocess.run([sys.executable, '-m', 'pip', 'install', '{}/.'.format(PathManager.get_package_path()),
                        '--upgrade'])

        cls.package_path = os.path.join(PathManager.get_package_path())
        cls.test_path = os.path.join(PathManager.get_test_path())
        cls.outdir_path = os.path.join(cls.test_path, 'outdir')
        # during development of the test, this prevents errors
        shutil.rmtree(cls.outdir_path, ignore_errors=True)
        pathlib.Path(cls.outdir_path).mkdir(parents=True, exist_ok=True)
        os.environ['VTAM_LOG_VERBOSITY'] = str(10)

        ############################################################################################
        #
        # Download fastq test dataset
        #
        ############################################################################################

        sorted_tar_path = os.path.join(cls.outdir_path, "sorted.tar.gz")
        if not os.path.isfile(sorted_tar_path):
            urllib.request.urlretrieve(sorted_tar_gz_url, sorted_tar_path)
        tar = tarfile.open(sorted_tar_path, "r:gz")
        tar.extractall(path=cls.outdir_path)
        tar.close()

        ############################################################################################
        #
        # Paths
        #
        ############################################################################################

        cls.asvtable_path = os.path.join(cls.outdir_path, "asvtable_default.tsv")

        cls.args = {}
        cls.args['readinfo'] = os.path.join(os.path.dirname(__file__), "readinfo.tsv")
        cls.args['params_lfn_variant'] = os.path.join(os.path.dirname(__file__), "params_lfn_variant.yml")
        cls.args['params_lfn_variant_replicate'] = os.path.join(os.path.dirname(__file__), "params_lfn_variant_replicate.yml")

    def test_filter_params_lfn_variant_replicate(self):

        ############################################################################################
        #
        # Wrong
        #
        ############################################################################################

        cmd = "vtam filter --db db.sqlite --readinfo {readinfo} --readdir sorted " \
              "--asvtable asvtable_default.tsv --params {params_lfn_variant} --until FilterLFN " \
              "--lfn_variant_replicate".format(**self.args)
        result = subprocess.run(shlex.split(cmd), cwd=self.outdir_path)

        self.assertEqual(result.returncode, 1)

        ############################################################################################
        #
        # Wrong
        #
        ############################################################################################

        cmd = "vtam filter --db db.sqlite --readinfo {readinfo} --readdir sorted " \
              "--asvtable asvtable_default.tsv --params {params_lfn_variant_replicate} --until FilterLFN " \
              "".format(**self.args)
        result = subprocess.run(shlex.split(cmd), cwd=self.outdir_path)

        self.assertEqual(result.returncode, 1)

        ############################################################################################
        #
        # Right
        #
        ############################################################################################

        cmd = "vtam filter --db db.sqlite --readinfo {readinfo} --readdir sorted " \
              "--asvtable asvtable_default.tsv --params {params_lfn_variant} --until FilterLFN " \
              "".format(**self.args)
        result = subprocess.run(shlex.split(cmd), cwd=self.outdir_path)

        self.assertEqual(result.returncode, 0)

        ############################################################################################
        #
        # Right
        #
        ############################################################################################

        cmd = "vtam filter --db db.sqlite --readinfo {readinfo} --readdir sorted " \
              "--asvtable asvtable_default.tsv --params {params_lfn_variant_replicate} --until FilterLFN " \
              "--lfn_variant_replicate".format(**self.args)
        result = subprocess.run(shlex.split(cmd), cwd=self.outdir_path)

        self.assertEqual(result.returncode, 0)

    @classmethod
    def tearDownClass(cls):

        shutil.rmtree(cls.outdir_path, ignore_errors=True)
