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
from vtam.utils import pip_install_vtam_for_tests
from vtam.utils.constants import fastq_tar_gz_url1, fastq_tar_gz_url2, fastq_tar_gz_url3
from tqdm import tqdm
from vtam.utils import tqdm_hook


@unittest.skipUnless(not sys.platform.startswith("win"), "Test does not work with Windows")
# Not working with windows because of commands in snake.tuto.data
class TestTutorialSnakemakeOneMarker(unittest.TestCase):

    """Will test snakemake with makeknownoccurrences (Castet)"""

    @classmethod
    def setUpClass(cls):

        ########################################################################
        #
        # These tests need the vtam command in the path
        #
        ########################################################################

        pip_install_vtam_for_tests()  # vtam needs to be in the path

        cls.package_path = PathManager.get_package_path()
        cls.test_path = PathManager.get_test_path()

        cls.outdir_path = os.path.join(cls.test_path, 'outdir')
        shutil.rmtree(cls.outdir_path, ignore_errors=True)
        cls.outdir_data_path = os.path.join(cls.outdir_path, 'data')
        pathlib.Path(cls.outdir_data_path).mkdir(parents=True, exist_ok=True)

        cls.outdir_download_path = os.path.join(cls.test_path, 'outdir_download')
        pathlib.Path(cls.outdir_download_path).mkdir(parents=True, exist_ok=True)

        cls.snakefile_tuto_data = os.path.join(cls.package_path, "data/snake.tuto.data_makeknownoccurrences.yml")

        ############################################################################################
        #
        # Set command args
        #
        ############################################################################################

        cls.args = {}
        cls.args['package_path'] = cls.package_path
        cls.args['snake_tuto_data'] = cls.snakefile_tuto_data

        ############################################################################################
        #
        # Download fastq test dataset
        #
        ############################################################################################

        fastq_tar_path = os.path.join(cls.outdir_download_path, "fastq.tar.gz")
        # Test first in local dir, otherwise in the remote URLs
        if not os.path.isfile(fastq_tar_path) or pathlib.Path(fastq_tar_path).stat().st_size < 1000000:
            try:
                # urllib.request.urlretrieve(fastq_tar_gz_url1, fastq_tar_path, MyProgressBar())
                with tqdm(...) as t:
                    t.set_description(os.path.basename(fastq_tar_path))
                    urllib.request.urlretrieve(fastq_tar_gz_url1, fastq_tar_path, reporthook=tqdm_hook(t))
            except Exception:
                try:
                    # urllib.request.urlretrieve(fastq_tar_gz_url2, fastq_tar_path, MyProgressBar())
                    with tqdm(...) as t:
                        t.set_description(os.path.basename(fastq_tar_path))
                        urllib.request.urlretrieve(fastq_tar_gz_url2, fastq_tar_path, reporthook=tqdm_hook(t))
                except Exception:
                    # urllib.request.urlretrieve(fastq_tar_gz_url3, fastq_tar_path, MyProgressBar())
                    with tqdm(...) as t:
                        t.set_description(os.path.basename(fastq_tar_path))
                        urllib.request.urlretrieve(fastq_tar_gz_url3, fastq_tar_path, reporthook=tqdm_hook(t))
        tar = tarfile.open(fastq_tar_path, "r:gz")
        tar.extractall(path=cls.outdir_path)
        tar.close()

        ############################################################################################
        #
        # Copy data to directory tree
        #
        ############################################################################################

        cmd = "snakemake --cores 1 -s {snake_tuto_data} --config MARKER=mfzr " \
              "PROJECT=asper1 PACKAGE_PATH={package_path} --until all_one_marker_makeknownoccurrences".format(**cls.args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=cls.outdir_path)

    def test_01_mfzr_makeknownoccurrences(self):

        snakeconfig = os.path.join("asper1", "user_input", "snakeconfig_mfzr_makeknownoccurrences.yml")
        cmd = "snakemake --printshellcmds --resources db=1 --snakefile snakefile_makeknownoccurrences.yml --cores 4 --configfile {} --until makeknownoccurrences".format(snakeconfig)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=self.outdir_path)

        known_occurrences_mfzr = os.path.join(self.outdir_path, "asper1", "run1_mfzr", "known_occurrences.tsv")
        known_occurrences_mfzr_bak = os.path.join(self.test_path, "test_files", "known_occurrences_mfzr.tsv")
        self.assertTrue(filecmp.cmp(known_occurrences_mfzr, known_occurrences_mfzr_bak, shallow=True))

    def test_02_mfzr_makeknownoccurrences_optimize(self):

        snakeconfig = os.path.join("asper1", "user_input", "snakeconfig_mfzr_makeknownoccurrences.yml")
        cmd = "snakemake --printshellcmds --resources db=1 --snakefile snakefile_makeknownoccurrences.yml --cores 4 --configfile {} --until optimize".format(snakeconfig)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=self.outdir_path)

    @classmethod
    def tearDownClass(cls):

        shutil.rmtree(cls.outdir_path, ignore_errors=True)
