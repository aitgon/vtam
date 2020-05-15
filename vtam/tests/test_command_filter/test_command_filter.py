import filecmp
import os
import pathlib
import shlex
import shutil
import subprocess
import sys
import tarfile
import unittest
import urllib

from vtam.utils.PathManager import PathManager
from vtam.utils.constants import sorted_tar_gz_url


class TestWrapperFilterLFNCommands(unittest.TestCase):

    """Will test main commands based on a complete test dataset"""

    def setUp(self):

        # vtam needs to be in the tsv_path
        subprocess.run([sys.executable, '-m', 'pip', 'install', '{}/.'.format(PathManager.get_package_path()),
                        '--upgrade'])

        self.test_outdir_path = os.path.join(PathManager.get_test_path(), 'outdir')
        shutil.rmtree(self.test_outdir_path, ignore_errors=True)  # during development of the test, this prevents errors
        pathlib.Path(self.test_outdir_path).mkdir(parents=True, exist_ok=True)

        ################################################################################################################
        #
        # Download fastq test dataset
        #
        ################################################################################################################

        sorted_tar_path = os.path.join(self.test_outdir_path, "sorted.tar.gz")
        if not os.path.isfile(sorted_tar_path):
            urllib.request.urlretrieve(sorted_tar_gz_url, sorted_tar_path)
        tar = tarfile.open(sorted_tar_path, "r:gz")
        tar.extractall(path=self.test_outdir_path)
        tar.close()

        ################################################################################################################
        #
        # Paths
        #
        ################################################################################################################

        self.sorted_dir_path = os.path.join(self.test_outdir_path, "sorted")
        self.sortedreadinfo_path = os.path.join(self.sorted_dir_path, "readinfo.tsv")

        self.log_path = os.path.join(self.test_outdir_path, "vtam.log")

        self.asvtable_path = os.path.join(self.test_outdir_path, "asvtable_default.tsv")

        self.args = {}
        self.args['sorted'] = self.sorted_dir_path
        self.args['db'] = os.path.join(self.test_outdir_path, "db.sqlite")
        self.args['readinfo'] = self.sortedreadinfo_path
        self.args['readdir'] = self.sorted_dir_path
        self.args['asvtable'] = self.asvtable_path
        self.args['log'] = self.log_path


    def test_01(self):

        ################################################################################################################
        #
        # Command Filter
        #
        ################################################################################################################

        cmd = "vtam filter --db {db} --readinfo {readinfo} --readdir {readdir} --asvtable {asvtable} " \
              "-v --log {log}".format(**self.args)
        subprocess.run(shlex.split(cmd))

        asvtable_bak_path = os.path.join(os.path.dirname(__file__), "asvtable_default.tsv")
        self.assertTrue(filecmp.cmp(self.asvtable_path, asvtable_bak_path, shallow=False))

    def tearDown(self):

        subprocess.run([sys.executable, '-m', 'pip', 'uninstall', 'vtam', '-y'])
        shutil.rmtree(self.test_outdir_path, ignore_errors=True)
