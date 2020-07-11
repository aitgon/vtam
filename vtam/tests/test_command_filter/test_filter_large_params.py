from urllib import request
from vtam.utils.PathManager import PathManager
from vtam.utils.constants import sorted_dryad_f40v5_tar_gz_url
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


@unittest.skipIf(request.urlopen(sorted_dryad_f40v5_tar_gz_url).getcode() != 200,
                 "This test requires an internet connection!")
class TestFilterLargeParams(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.package_path = os.path.join(PathManager.get_package_path())
        # vtam needs to be in the tsv_path
        subprocess.run([sys.executable, '-m', 'pip', 'install', '.', '--upgrade'], cwd=cls.package_path)

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

        sorted_tar_path = os.path.join(cls.outdir_path, "sorted.tar.gz")
        if not os.path.isfile(sorted_tar_path):
            urllib.request.urlretrieve(sorted_dryad_f40v5_tar_gz_url, sorted_tar_path)
        tar = tarfile.open(sorted_tar_path, "r:gz")
        tar.extractall(path=cls.outdir_data_path)
        tar.close()

        ############################################################################################
        #
        # Args
        #
        ############################################################################################

        cls.args = {}
        cls.args['params'] = os.path.join(os.path.dirname(__file__), "params_large.yml")
        cls.args['asvtable_default'] = os.path.join(cls.outdir_data_path, "asvtable_default.tsv")
        cls.args['blastdbdir'] = os.path.join(cls.outdir_data_path, "coi_blast_db_dir")
        cls.args['taxonomy'] = os.path.join(cls.outdir_data_path, "taxonomy.tsv")

        ############################################################################################
        #
        # Run 'vtam taxonomy'
        #
        ############################################################################################

        cmd = "vtam taxonomy --output {taxonomy} --precomputed".format(**cls.args)
        # print(cmd)
        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=cls.outdir_data_path)

        ############################################################################################
        #
        # Run 'vtam coi_blast_db'
        #
        ############################################################################################

        cmd = "vtam coi_blast_db --blastdbdir {blastdbdir} --blastdbname coi_blast_db_20191211".format(**cls.args)
        # print(cmd)
        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=cls.outdir_data_path)

        ############################################################################################
        #
        # Command Filter
        #
        ############################################################################################

        cmd = "vtam filter --lfn_variant_replicate --db db.sqlite --readinfo sorted/readinfo.tsv " \
              "--readdir sorted --asvtable {asvtable_default} --params {params}".format(**cls.args)
        # print(cmd)
        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, cwd=cls.outdir_data_path)
        pass

    def setUp(self):

        self.outdir_thistest_path = os.path.join(self.outdir_path, 'thistest')
        # during development of the test, this prevents errors
        pathlib.Path(self.outdir_thistest_path).mkdir(parents=True, exist_ok=True)
        os.environ['VTAM_LOG_VERBOSITY'] = str(10)

    def test_01_filter(self):

        asvtable_path = self.args['asvtable_default']
        asvtable_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5/asvtable_default.tsv")
        self.assertTrue(filecmp.cmp(asvtable_path, asvtable_bak_path, shallow=False))

    def test_02_taxasign(self):

        ############################################################################################
        #
        # Command Filter
        # #
        ############################################################################################

        cmd = "vtam taxassign --variants {asvtable_default} --output asvtable_default_taxa.tsv " \
                  "--db db.sqlite --blastdbdir {blastdbdir} --blastdbname coi_blast_db_20191211 --taxonomy {taxonomy}".format(**self.args)
        # print(cmd)
        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, cwd=self.outdir_thistest_path)

        asvtable_path = os.path.join(self.outdir_thistest_path, "asvtable_default_taxa.tsv")
        asvtable_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5/asvtable_default_taxa.tsv")
        self.assertTrue(filecmp.cmp(asvtable_path, asvtable_bak_path, shallow=False))

    def tearDown(self):

        shutil.rmtree(self.outdir_thistest_path, ignore_errors=True)

    @classmethod
    def tearDownClass(cls):

        shutil.rmtree(cls.outdir_path, ignore_errors=True)
