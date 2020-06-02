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

        cls.sorted_dir_path = os.path.join(cls.outdir_path, "sorted")
        cls.sortedreadinfo_path = os.path.join(cls.sorted_dir_path, "readinfo.tsv")
        cls.asvtable_path = os.path.join(cls.outdir_path, "asvtable_default.tsv")

        cls.known_occurrences = os.path.join(cls.package_path, "doc/data/known_occurrences.tsv")

        cls.log_path = os.path.join(cls.outdir_path, "vtam.log")

        cls.args = {}
        cls.args['db'] = os.path.join(cls.outdir_path, "db.sqlite")
        cls.args['readinfo'] = cls.sortedreadinfo_path
        cls.args['readdir'] = cls.sorted_dir_path
        cls.args['known_occurrences'] = cls.known_occurrences
        cls.args['outdir'] = cls.outdir_path
        cls.args['log'] = cls.log_path
        cls.args['asvtable'] = cls.asvtable_path

        ############################################################################################
        #
        # Command Optimize
        #
        ############################################################################################

        cmd = "vtam filter --db {db} --readinfo {readinfo} --readdir {readdir} " \
              "--asvtable {asvtable} -v --log {log} --until VariantReadCount".format(**cls.args)
        subprocess.run(shlex.split(cmd))

    def test_01_filter_lfn_variant(self):

        ################################################################################################################
        #
        # Command Filter
        #
        ################################################################################################################

        cmd = "vtam filter --db {db} --readinfo {readinfo} --readdir {readdir} --asvtable {asvtable} " \
              "-v --log {log}".format(**self.args)
        subprocess.run(shlex.split(cmd))

        asvtable_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/asvtable_default.tsv")
        self.assertTrue(filecmp.cmp(self.asvtable_path, asvtable_bak_path, shallow=False))

    def test_01_filter_lfn_variant_replicate(self):

        ################################################################################################################
        #
        # Command Filter
        #
        ################################################################################################################

        cmd = "vtam filter --db {db} --readinfo {readinfo} --readdir {readdir} --asvtable {asvtable} --lfn_variant_replicate " \
              "-v --log {log}".format(**self.args)
        subprocess.run(shlex.split(cmd))

        asvtable_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/asvtable_default_lfn_variant_replicate.tsv")
        self.assertTrue(filecmp.cmp(self.asvtable_path, asvtable_bak_path, shallow=False))

    def test_02_optimize_lfn_biosample_replicate(self):

        cmd = "vtam optimize --db {db} --readinfo {readinfo} --readdir {readdir} " \
              "--known_occurrences {known_occurrences} --outdir {outdir} --until OptimizeLFNbiosampleReplicate " \
              "-v --log {log}".format(**self.args)
        subprocess.run(shlex.split(cmd))

        optimize_lfn_biosample_replicate_path = os.path.join(self.outdir_path, "optimize_lfn_biosample_replicate.tsv")
        optimize_lfn_biosample_replicate_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/optimize_lfn_biosample_replicate.tsv")
        self.assertTrue(filecmp.cmp(optimize_lfn_biosample_replicate_path, optimize_lfn_biosample_replicate_bak_path, shallow=False))

    def test_03_optimize_pcr_error(self):

        cmd = "vtam optimize --db {db} --readinfo {readinfo} --readdir {readdir} " \
              "--known_occurrences {known_occurrences} --outdir {outdir} --until OptimizePCRerror " \
              "-v --log {log}".format(**self.args)
        subprocess.run(shlex.split(cmd))

        optimize_pcr_error_path = os.path.join(self.outdir_path, "optimize_pcr_error.tsv")
        optimize_pcr_error_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/optimize_pcr_error.tsv")
        self.assertTrue(filecmp.cmp(optimize_pcr_error_path, optimize_pcr_error_bak_path, shallow=False))

    def test_04_optimize_lfn_read_count_variant(self):

        cmd = "vtam optimize --db {db} --readinfo {readinfo} --readdir {readdir} " \
              "--known_occurrences {known_occurrences} --outdir {outdir} --until OptimizeLFNreadCountAndLFNvariant " \
              "-v --log {log}".format(**self.args)
        subprocess.run(shlex.split(cmd))

        optimize_lfn_read_count_variant_path = os.path.join(self.outdir_path,
                                                                "optimize_lfn_read_count_and_lfn_variant.tsv")
        optimize_lfn_read_count_variant_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/optimize_lfn_read_count_and_lfn_variant.tsv")
        self.assertTrue(filecmp.cmp(optimize_lfn_read_count_variant_path, optimize_lfn_read_count_variant_bak_path, shallow=False))

        optimize_lfn_variant_specific_path = os.path.join(self.outdir_path,
                                                                "optimize_lfn_variant_specific.tsv")
        optimize_lfn_variant_specific_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/optimize_lfn_variant_specific.tsv")
        self.assertTrue(filecmp.cmp(optimize_lfn_variant_specific_path, optimize_lfn_variant_specific_bak_path, shallow=False))

    def test_05_optimize_lfn_read_count_variant_replicate(self):

        cmd = "vtam optimize --db {db} --readinfo {readinfo} --readdir {readdir} " \
              "--known_occurrences {known_occurrences} --outdir {outdir} --until OptimizeLFNreadCountAndLFNvariant --lfn_variant_replicate " \
              "-v --log {log}".format(**self.args)
        subprocess.run(shlex.split(cmd))

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
