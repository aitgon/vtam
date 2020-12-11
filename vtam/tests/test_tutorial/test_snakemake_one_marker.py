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

    """Will test main commands based on a complete test dataset"""

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

        cls.snakefile_tuto_data = os.path.join(cls.package_path, "data/snake.tuto.data.yml")

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
              "PROJECT=asper1 PACKAGE_PATH={package_path} --until all_one_marker".format(**cls.args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=cls.outdir_path)

        cmd = "snakemake --cores 1 -s {snake_tuto_data} --config MARKER=zfzr " \
                  "PROJECT=asper1 PACKAGE_PATH={package_path} --until all_one_marker".format(**cls.args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=cls.outdir_path)

    def test_01_mfzr_filter(self):

        snakeconfig = os.path.join("asper1", "user_input", "snakeconfig_mfzr.yml")
        cmd = "snakemake --printshellcmds --resources db=1 --snakefile snakefile.yml --cores 4 --configfile {} --until asvtable_taxa".format(snakeconfig)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=self.outdir_path)

    def test_02_mfzr_optimize(self):

        snakeconfig = os.path.join("asper1", "user_input", "snakeconfig_mfzr.yml")
        cmd = "snakemake --printshellcmds --resources db=1 --snakefile snakefile.yml --cores 4 --configfile {} --until optimize".format(snakeconfig)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=self.outdir_path)

    def test_03_mfzr_filter_optimized(self):

        snakeconfig = os.path.join("asper1", "user_input", "snakeconfig_mfzr.yml")
        cmd = "snakemake --printshellcmds --resources db=1 --snakefile snakefile.yml --cores 4 --configfile {} --until asvtable_optimized_taxa".format(snakeconfig)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=self.outdir_path)

    def test_04_zfzr_filter(self):

        snakeconfig = os.path.join("asper1", "user_input", "snakeconfig_zfzr.yml")
        cmd = "snakemake --printshellcmds --resources db=1 --snakefile snakefile.yml --cores 4 --configfile {} --until asvtable_taxa".format(snakeconfig)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=self.outdir_path)

    def test_05_zfzr_optimize(self):

        snakeconfig = os.path.join("asper1", "user_input", "snakeconfig_zfzr.yml")
        cmd = "snakemake --printshellcmds --resources db=1 --snakefile snakefile.yml --cores 4 --configfile {} --until optimize".format(snakeconfig)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=self.outdir_path)

    def test_06_zfzr_filter_optimized(self):

        snakeconfig = os.path.join("asper1", "user_input", "snakeconfig_zfzr.yml")
        cmd = "snakemake --printshellcmds --resources db=1 --snakefile snakefile.yml --cores 4 --configfile {} --until asvtable_optimized_taxa".format(snakeconfig)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=self.outdir_path)

    def test_07_pool(self):

        db = os.path.join("asper1", "db.sqlite")
        runmarker = os.path.join("asper1", "user_input", "pool_run_marker.tsv")
        asvtable_pooled = os.path.join("asper1", "asvtable_pooled_mfzr_zfzr.tsv")
        log = os.path.join("asper1", "vtam.log")
        args = {'db': db, 'runmarker': runmarker, 'asvtable_pooled': asvtable_pooled, 'log': log}
        cmd = "vtam pool --db asper1/db.sqlite --runmarker {runmarker} --asvtable {asvtable_pooled} --log {log} -v".format(**args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=self.outdir_path)

    def test_08_pool_taxa(self):

        db = os.path.join("asper1", "db.sqlite")
        asvtable = os.path.join("asper1", "asvtable_pooled_mfzr_zfzr.tsv")
        output = os.path.join("asper1", "asvtable_pooled_mfzr_zfzr_taxa.tsv")
        taxonomy = os.path.join("vtam_db", "taxonomy.tsv")
        blastdbdir = os.path.join("vtam_db", "coi_blast_db")
        log = os.path.join("asper1", "vtam.log")
        args = {'db': db, 'asvtable': asvtable, 'output': output, 'taxonomy': taxonomy, 'blastdbdir': blastdbdir, 'log': log}
        cmd = "vtam taxassign --db {db} --asvtable {asvtable} " \
              "--output {output} --taxonomy {taxonomy} " \
              "--blastdbdir {blastdbdir} --blastdbname coi_blast_db_20200420 --log {log} -v".format(**args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=self.outdir_path)

    @classmethod
    def tearDownClass(cls):

        shutil.rmtree(cls.outdir_path, ignore_errors=True)
