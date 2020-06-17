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
        # shutil.rmtree(cls.outdir_path, ignore_errors=True)
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
        cls.args['known_occurrences'] = os.path.join(cls.package_path, "doc/data/known_occurrences.tsv")
        cls.args['lfn_variant_cutoff_specific'] = os.path.join(
            cls.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/optimize_lfn_variant_specific.tsv")
        cls.args['lfn_variant_replicate_cutoff_specific'] = os.path.join(
            cls.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/optimize_lfn_variant_replicate_specific.tsv")

        ############################################################################################
        #
        # Command Optimize
        #
        ############################################################################################

        cmd = "vtam filter --db db.sqlite --readinfo sorted/readinfo.tsv --readdir sorted " \
              "--asvtable asvtable_default.tsv -v --log vtam.log --until VariantReadCount".format(**cls.args)
        subprocess.run(shlex.split(cmd), cwd=cls.outdir_path)

    def test_01_filter_lfn_variant(self):

        cmd = "vtam filter --db db.sqlite --readinfo sorted/readinfo.tsv --readdir sorted --asvtable asvtable_default.tsv " \
              "-v --log vtam.log".format(**self.args)
        subprocess.run(shlex.split(cmd), cwd=self.outdir_path)

        asvtable_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/asvtable_default.tsv")
        self.assertTrue(filecmp.cmp(self.asvtable_path, asvtable_bak_path, shallow=False))

    def test_02_filter_lfn_variant_replicate(self):

        cmd = "vtam filter --db db.sqlite --readinfo sorted/readinfo.tsv --readdir sorted --asvtable asvtable_default.tsv " \
              "--lfn_variant_replicate -v --log vtam.log".format(**self.args)
        subprocess.run(shlex.split(cmd), cwd=self.outdir_path)

        asvtable_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/asvtable_default_lfn_variant_replicate.tsv")
        self.assertTrue(filecmp.cmp(self.asvtable_path, asvtable_bak_path, shallow=False))

    def test_03_filter_lfn_variant_cutoff_specific(self):

        cmd = "vtam filter --db db.sqlite --readinfo sorted/readinfo.tsv --readdir sorted --asvtable asvtable_default.tsv " \
              "-v --log vtam.log --cutoff_specific {lfn_variant_cutoff_specific}".format(**self.args)
        subprocess.run(shlex.split(cmd), cwd=self.outdir_path)

        asvtable_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/asvtable_default_lfn_variant_cutoff_specific.tsv")
        self.assertTrue(filecmp.cmp(self.asvtable_path, asvtable_bak_path, shallow=False))

    def test_04_filter_lfn_variant_replicate_cutoff_specific(self):

        cmd = "vtam filter --db db.sqlite --readinfo sorted/readinfo.tsv --readdir sorted --asvtable asvtable_default.tsv " \
              "--lfn_variant_replicate -v --log vtam.log --cutoff_specific {lfn_variant_replicate_cutoff_specific}".format(**self.args)
        subprocess.run(shlex.split(cmd), cwd=self.outdir_path)

        asvtable_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/asvtable_default_lfn_variant_replicate_cutoff_specific.tsv")
        self.assertTrue(filecmp.cmp(self.asvtable_path, asvtable_bak_path, shallow=False))

    def test_05_optimize_lfn_biosample_replicate(self):

        cmd = "vtam optimize --db db.sqlite --readinfo sorted/readinfo.tsv --readdir sorted " \
              "--known_occurrences {known_occurrences} --outdir . --until OptimizeLFNbiosampleReplicate " \
              "-v --log vtam.log".format(**self.args)
        subprocess.run(shlex.split(cmd), cwd=self.outdir_path)

        optimize_lfn_biosample_replicate_path = os.path.join(self.outdir_path, "optimize_lfn_biosample_replicate.tsv")
        optimize_lfn_biosample_replicate_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/optimize_lfn_biosample_replicate.tsv")
        self.assertTrue(filecmp.cmp(optimize_lfn_biosample_replicate_path, optimize_lfn_biosample_replicate_bak_path, shallow=False))

    def test_06_optimize_pcr_error(self):

        cmd = "vtam optimize --db db.sqlite --readinfo sorted/readinfo.tsv --readdir sorted " \
              "--known_occurrences {known_occurrences} --outdir . --until OptimizePCRerror " \
              "-v --log vtam.log".format(**self.args)
        subprocess.run(shlex.split(cmd), cwd=self.outdir_path)

        optimize_pcr_error_path = os.path.join(self.outdir_path, "optimize_pcr_error.tsv")
        optimize_pcr_error_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/optimize_pcr_error.tsv")
        self.assertTrue(filecmp.cmp(optimize_pcr_error_path, optimize_pcr_error_bak_path, shallow=False))

    def test_07_optimize_lfn_read_count_variant(self):

        cmd = "vtam optimize --db db.sqlite --readinfo sorted/readinfo.tsv --readdir sorted " \
              "--known_occurrences {known_occurrences} --outdir . --until OptimizeLFNreadCountAndLFNvariant " \
              "-v --log vtam.log".format(**self.args)

        subprocess.run(shlex.split(cmd), cwd=self.outdir_path)
        optimize_lfn_read_count_variant_path = os.path.join(self.outdir_path,
                                                                "optimize_lfn_read_count_and_lfn_variant.tsv")
        optimize_lfn_read_count_variant_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/optimize_lfn_read_count_and_lfn_variant.tsv")
        self.assertTrue(filecmp.cmp(optimize_lfn_read_count_variant_path, optimize_lfn_read_count_variant_bak_path, shallow=False))

        optimize_lfn_variant_specific_path = os.path.join(self.outdir_path,
                                                                "optimize_lfn_variant_specific.tsv")
        optimize_lfn_variant_specific_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/optimize_lfn_variant_specific.tsv")
        self.assertTrue(filecmp.cmp(optimize_lfn_variant_specific_path, optimize_lfn_variant_specific_bak_path, shallow=False))

    def test_08_optimize_lfn_read_count_variant_replicate(self):

        cmd = "vtam optimize --db db.sqlite --readinfo sorted/readinfo.tsv --readdir sorted " \
              "--known_occurrences {known_occurrences} --outdir . --until OptimizeLFNreadCountAndLFNvariant --lfn_variant_replicate " \
              "-v --log vtam.log".format(**self.args)
        subprocess.run(shlex.split(cmd), cwd=self.outdir_path)
        optimize_lfn_read_count_variant_replicate_path = os.path.join(self.outdir_path,
                                                                "optimize_lfn_read_count_and_lfn_variant_replicate.tsv")
        optimize_lfn_read_count_variant_replicate_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/optimize_lfn_read_count_and_lfn_variant_replicate.tsv")
        self.assertTrue(filecmp.cmp(optimize_lfn_read_count_variant_replicate_path, optimize_lfn_read_count_variant_replicate_bak_path, shallow=False))

        optimize_lfn_variant_replicate_specific_path = os.path.join(self.outdir_path,
                                                                "optimize_lfn_variant_replicate_specific.tsv")
        optimize_lfn_variant_replicate_specific_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/optimize_lfn_variant_replicate_specific.tsv")
        self.assertTrue(filecmp.cmp(optimize_lfn_variant_replicate_specific_path, optimize_lfn_variant_replicate_specific_bak_path, shallow=False))


    @classmethod
    def tearDownClass(cls):

        shutil.rmtree(cls.outdir_path, ignore_errors=True)
