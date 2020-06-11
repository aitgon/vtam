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
        cls.args['params'] = os.path.join(os.path.dirname(__file__), "params.yml")

    def test_01_filter_params_no_params_consecutively(self):

        ############################################################################################
        #
        # One replicate, variants left, because min_replicate_nb =1
        #
        ############################################################################################

        cmd = "vtam filter --db db.sqlite --readinfo {readinfo} --readdir sorted " \
              "--asvtable asvtable_default.tsv --params {params} --until FilterMinReplicateNumber".format(**self.args)
        subprocess.run(shlex.split(cmd), cwd=self.outdir_path)

        db_path = os.path.join(self.outdir_path, "db.sqlite")
        con = sqlite3.connect(db_path)
        cur = con.cursor()
        cur_result = cur.execute('SELECT COUNT(*) from FilterMinReplicateNumber where filter_delete=0').fetchone()
        self.assertTrue(cur_result[0] > 0)
        con.close()

        ############################################################################################
        #
        # One replicate, no variants left, because default min_replicate_nb = 2
        #
        ############################################################################################

        cmd = "vtam filter --db db.sqlite --readinfo {readinfo} --readdir sorted " \
              "--asvtable asvtable_default.tsv --until FilterMinReplicateNumber".format(**self.args)
        subprocess.run(shlex.split(cmd), cwd=self.outdir_path)

        db_path = os.path.join(self.outdir_path, "db.sqlite")
        con = sqlite3.connect(db_path)
        cur = con.cursor()
        cur_result = cur.execute('SELECT COUNT(*) from FilterMinReplicateNumber where filter_delete=0').fetchone()
        self.assertFalse(cur_result[0] > 0)
        con.close()

        ############################################################################################
        #
        # One replicate, variants left, because min_replicate_nb =1
        #
        ############################################################################################

        cmd = "vtam filter --db db.sqlite --readinfo {readinfo} --readdir sorted " \
              "--asvtable asvtable_default.tsv --params {params} --until FilterMinReplicateNumber".format(**self.args)
        subprocess.run(shlex.split(cmd), cwd=self.outdir_path)

        db_path = os.path.join(self.outdir_path, "db.sqlite")
        con = sqlite3.connect(db_path)
        cur = con.cursor()
        cur_result = cur.execute('SELECT COUNT(*) from FilterMinReplicateNumber where filter_delete=0').fetchone()
        self.assertTrue(cur_result[0] > 0)
        con.close()

    @classmethod
    def tearDownClass(cls):

        shutil.rmtree(cls.outdir_path, ignore_errors=True)
