import filecmp
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
from tqdm import tqdm
from vtam.utils import tqdm_hook


class TestCommandsFilter(unittest.TestCase):

    """Will test main commands based on a complete test dataset"""

    def setUp(self):

        self.package_path = os.path.join(PathManager.get_package_path())
        self.test_path = os.path.join(PathManager.get_test_path())
        self.outdir_path = os.path.join(self.test_path, 'outdir')
        self.outdir_data_path = os.path.join(self.outdir_path, 'data')
        shutil.rmtree(self.outdir_path, ignore_errors=True)
        pathlib.Path(self.outdir_data_path).mkdir(parents=True, exist_ok=True)

        ############################################################################################
        #
        # Download test dataset
        #
        ############################################################################################

        sorted_tar_path = os.path.join(self.outdir_data_path, "sorted.tar.gz")
        # Test first in local dir, otherwise in the remote URLs
        if not os.path.isfile(sorted_tar_path) or pathlib.Path(sorted_tar_path).stat().st_size < 1000000:
            try:
                # urllib.request.urlretrieve(sorted_tar_gz_url1, sorted_tar_path, MyProgressBar())
                with tqdm(...) as t:
                    t.set_description(os.path.basename(sorted_tar_path))
                    urllib.request.urlretrieve(sorted_tar_gz_url1, sorted_tar_path, reporthook=tqdm_hook(t))
            except Exception:
                try:
                    # urllib.request.urlretrieve(sorted_tar_gz_url2, sorted_tar_path, MyProgressBar())
                    with tqdm(...) as t:
                        t.set_description(os.path.basename(sorted_tar_path))
                        urllib.request.urlretrieve(sorted_tar_gz_url2, sorted_tar_path, reporthook=tqdm_hook(t))
                except Exception:
                    # urllib.request.urlretrieve(sorted_tar_gz_url3, sorted_tar_path, MyProgressBar())
                    with tqdm(...) as t:
                        t.set_description(os.path.basename(sorted_tar_path))
                        urllib.request.urlretrieve(sorted_tar_gz_url3, sorted_tar_path, reporthook=tqdm_hook(t))
        tar = tarfile.open(sorted_tar_path, "r:gz")
        tar.extractall(path=self.outdir_data_path)
        tar.close()

        self.outdir_thistest_path = os.path.join(self.outdir_path, 'thistest')
        # during development of the test, this prevents errors
        pathlib.Path(self.outdir_thistest_path).mkdir(parents=True, exist_ok=True)
        os.environ['VTAM_LOG_VERBOSITY'] = str(10)

        ############################################################################################
        #
        # Paths
        #
        ############################################################################################

        self.asvtable_path = os.path.join(self.outdir_thistest_path, "asvtable_default.tsv")

        self.args = {}
        self.args['sorteddir'] = os.path.join(self.outdir_data_path, 'sorted')
        self.args['sortedinfo'] = os.path.join(self.args['sorteddir'], "sortedinfo.tsv")
        self.args['known_occurrences'] = os.path.join(self.package_path, "data/example/known_occurrences.tsv")
        self.args['lfn_variant_cutoff_specific'] = os.path.join(
            self.test_path, "test_files_dryad.f40v5_small", "run1_mfzr_zfzr", "optimize_lfn_variant_specific.tsv")
        self.args['lfn_variant_replicate_cutoff_specific'] = os.path.join(
            self.test_path, "test_files_dryad.f40v5_small", "run1_mfzr_zfzr", "optimize_lfn_variant_replicate_specific.tsv")
        self.args['asvtable_default'] = self.asvtable_path

    def test_01_filter_lfn_variant(self):

        cmd = "vtam filter --db db.sqlite --sortedinfo {sortedinfo} --sorteddir {sorteddir} --asvtable {asvtable_default} " \
              "--log vtam.log".format(**self.args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, cwd=self.outdir_thistest_path)

        asvtable_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5_small", "run1_mfzr_zfzr", "asvtable_default.tsv")
        self.assertTrue(filecmp.cmp(self.asvtable_path, asvtable_bak_path, shallow=False))

    def test_02_filter_lfn_variant_replicate(self):

        cmd = "vtam filter --db db.sqlite --sortedinfo {sortedinfo} --sorteddir {sorteddir} --asvtable {asvtable_default} " \
              "--lfn_variant_replicate --log vtam.log".format(**self.args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, cwd=self.outdir_thistest_path)

        asvtable_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5_small", "run1_mfzr_zfzr", "asvtable_default_lfn_variant_replicate.tsv")
        self.assertTrue(filecmp.cmp(self.asvtable_path, asvtable_bak_path, shallow=False))

    def test_03_filter_lfn_variant_cutoff_specific(self):

        cmd = "vtam filter --db db.sqlite --sortedinfo {sortedinfo} --sorteddir {sorteddir} --asvtable {asvtable_default} " \
              "--log vtam.log --cutoff_specific {lfn_variant_cutoff_specific}".format(**self.args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, cwd=self.outdir_thistest_path)

        asvtable_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5_small", "run1_mfzr_zfzr", "asvtable_default_lfn_variant_cutoff_specific.tsv")
        self.assertTrue(filecmp.cmp(self.asvtable_path, asvtable_bak_path, shallow=False))

    def test_04_filter_lfn_variant_replicate_cutoff_specific(self):

        cmd = "vtam filter --db db.sqlite --sortedinfo {sortedinfo} --sorteddir {sorteddir} --asvtable {asvtable_default} " \
              "--lfn_variant_replicate --log vtam.log --cutoff_specific {lfn_variant_replicate_cutoff_specific}".format(**self.args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, cwd=self.outdir_thistest_path)

        asvtable_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5_small", "run1_mfzr_zfzr", "asvtable_default_lfn_variant_replicate_cutoff_specific.tsv")
        self.assertTrue(filecmp.cmp(self.asvtable_path, asvtable_bak_path, shallow=False))

    def test_09_filter_lfn_variant_with_known_occurrences(self):

        cmd = "vtam filter --db db.sqlite --sortedinfo {sortedinfo} --sorteddir {sorteddir} --asvtable {asvtable_default} " \
              "--known_occurrences {known_occurrences} " \
              "--log vtam.log".format(**self.args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, cwd=self.outdir_thistest_path)

        asvtable_bak_path = os.path.join(self.test_path, "test_files_dryad.f40v5_small", "run1_mfzr_zfzr", "asvtable_default_keep_occurrences.tsv")
        self.assertTrue(filecmp.cmp(self.asvtable_path, asvtable_bak_path, shallow=False))

    def tearDown(self):

        shutil.rmtree(self.outdir_thistest_path, ignore_errors=True)

        shutil.rmtree(self.outdir_path, ignore_errors=True)
