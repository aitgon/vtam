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
class TestFilterMinReplicateNumber(unittest.TestCase):

    """Will test main commands based on a complete test dataset"""

    def setUp(self):

        # vtam needs to be in the tsv_path
        cmd = '{} -m pip install . -q --upgrade'.format(sys.executable)
        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=PathManager.get_package_path())

        self.package_path = os.path.join(PathManager.get_package_path())
        self.test_path = os.path.join(PathManager.get_test_path())
        self.outdir_path = os.path.join(self.test_path, 'outdir')
        # during development of the test, this prevents errors
        shutil.rmtree(self.outdir_path, ignore_errors=True)
        pathlib.Path(self.outdir_path).mkdir(parents=True, exist_ok=True)
        os.environ['VTAM_LOG_VERBOSITY'] = str(10)

        ############################################################################################
        #
        # Download sorted fasta test dataset
        #
        ############################################################################################

        sorted_tar_path = os.path.join(self.outdir_path, "sorted.tar.gz")
        if not os.path.isfile(sorted_tar_path):
            urllib.request.urlretrieve(sorted_tar_gz_url, sorted_tar_path)
        tar = tarfile.open(sorted_tar_path, "r:gz")
        tar.extractall(path=self.outdir_path)
        tar.close()

        ############################################################################################
        #
        # Paths
        #
        ############################################################################################

        self.asvtable_path = os.path.join(self.outdir_path, "asvtable_default.tsv")

        self.args = {}
        self.args['readinfo'] = os.path.join(os.path.dirname(__file__), "readinfo.tsv")
        self.args['params'] = os.path.join(os.path.dirname(__file__), "params_min_replicate_number1.yml")
        self.args['params_lfn_variant'] = os.path.join(os.path.dirname(__file__), "params_lfn_variant.yml")
        self.args['params_lfn_variant_replicate'] = os.path.join(os.path.dirname(__file__), "params_lfn_variant_replicate.yml")

    def test_filter_params_no_params_consecutively(self):

        ############################################################################################
        #
        # One replicate, variants left, because min_replicate_nb =1
        #
        ############################################################################################

        cmd = "vtam filter --db db.sqlite --readinfo {readinfo} --readdir sorted " \
              "--asvtable asvtable_default.tsv --params {params} --until FilterMinReplicateNumber " \
              "--lfn_variant_replicate".format(**self.args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, cwd=self.outdir_path)

        db_path = os.path.join(self.outdir_path, "db.sqlite")
        con = sqlite3.connect(db_path)
        cur = con.cursor()
        cur_result = cur.execute('SELECT COUNT(*) from FilterMinReplicateNumber where filter_delete=0').fetchone()
        self.assertGreater(cur_result[0], 0)
        con.close()

        ############################################################################################
        #
        # One replicate, no variants left, because default min_replicate_nb = 2
        #
        ############################################################################################

        cmd = "vtam filter --db db.sqlite --readinfo {readinfo} --readdir sorted " \
              "--asvtable asvtable_default.tsv --until FilterMinReplicateNumber".format(**self.args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, cwd=self.outdir_path)

        db_path = os.path.join(self.outdir_path, "db.sqlite")
        con = sqlite3.connect(db_path)
        cur = con.cursor()
        cur_result = cur.execute('SELECT COUNT(*) from FilterMinReplicateNumber where filter_delete=0').fetchone()
        self.assertEqual(cur_result[0], 0)
        con.close()

        ############################################################################################
        #
        # One replicate, variants left, because min_replicate_nb =1
        #
        ############################################################################################

        cmd = "vtam filter --db db.sqlite --readinfo {readinfo} --readdir sorted " \
              "--asvtable asvtable_default.tsv --params {params} --until FilterMinReplicateNumber".format(**self.args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, cwd=self.outdir_path)

        db_path = os.path.join(self.outdir_path, "db.sqlite")
        con = sqlite3.connect(db_path)
        cur = con.cursor()
        cur_result = cur.execute('SELECT COUNT(*) from FilterMinReplicateNumber where filter_delete=0').fetchone()
        self.assertGreater(cur_result[0], 0)
        con.close()

    def tearDown(self):

        shutil.rmtree(self.outdir_path, ignore_errors=True)
