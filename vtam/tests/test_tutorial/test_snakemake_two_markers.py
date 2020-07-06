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

from vtam.utils.constants import fastq_tar_gz_url
from vtam.utils.PathManager import PathManager
from urllib import request


@unittest.skipIf(request.urlopen(fastq_tar_gz_url).getcode() != 200,
                 "This test requires an internet connection!")
class TestTutorialSnakemake(unittest.TestCase):

    """Will test main commands based on a complete test dataset"""

    @classmethod
    def setUpClass(cls):

        # vtam needs to be in the tsv_path
        subprocess.run([sys.executable, '-m', 'pip', 'install', '{}/.'.format(PathManager.get_package_path()),
                        '--upgrade'])

        cls.package_path = PathManager.get_package_path()
        cls.test_path = PathManager.get_test_path()
        cls.outdir_path = os.path.join(cls.test_path, 'outdir')
        # shutil.rmtree(cls.outdir_path, ignore_errors=True)
        pathlib.Path(cls.outdir_path).mkdir(parents=True, exist_ok=True)

        cls.snakefile_tuto_data = os.path.join(cls.package_path, "tools/snake.tuto.data.yml")

        ############################################################################################
        #
        # Set command args
        #
        ############################################################################################

        cls.args = {}
        cls.args['package_path'] = cls.package_path
        cls.args['snake_tuto_data'] = os.path.join(cls.package_path, "tools/snake.tuto.data.yml")

        ############################################################################################
        #
        # Download fastq test dataset
        #
        ############################################################################################

        fastq_tar_path = os.path.join(cls.outdir_path, "fastq.tar.gz")
        if not os.path.isfile(fastq_tar_path):
            urllib.request.urlretrieve(fastq_tar_gz_url, fastq_tar_path)
        tar = tarfile.open(fastq_tar_path, "r:gz")
        tar.extractall(path=cls.outdir_path)
        tar.close()

        ############################################################################################
        #
        # Copy data to directory tree
        #
        ############################################################################################

        cmd = "snakemake --cores 1 -s {snake_tuto_data} --config " \
                  "PROJECT=asper2 PACKAGE_PATH={package_path} --until all_two_markers".format(**cls.args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=cls.outdir_path)

    def test_01_filter(self):

        cmd = "snakemake --printshellcmds --resources db=1 --snakefile snakefile.yml --cores 4 --configfile asper2/user_input/snakeconfig.yml --until asvtable_taxa"

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=self.outdir_path)

    def test_02_optimize(self):

        cmd = "snakemake --printshellcmds --resources db=1 --snakefile snakefile.yml --cores 4 --configfile asper2/user_input/snakeconfig.yml --until optimize"

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=self.outdir_path)

        optimize_pcr_error = os.path.join(self.outdir_path, "asper2/run1/optimize_pcr_error.tsv")
        optimize_pcr_error_bak = os.path.join(self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/optimize_pcr_error.tsv")
        self.assertTrue(filecmp.cmp(optimize_pcr_error, optimize_pcr_error_bak, shallow=False))

    def test_03_mfzr_filter_optimized(self):

        cmd = "vtam filter --db asper2/db.sqlite --readinfo asper2/user_input/readinfo_mfzr.tsv --readdir asper2/run1/sorted --params asper2/user_input/params_mfzr.yml --asvtable asper2/run1/asvtable_params_mfzr.tsv -v --log asper2/vtam.log"

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=self.outdir_path)

    def test_03_mfzr_taxassign_optimized(self):

        cmd = "vtam taxassign --db asper2/db.sqlite --variants asper2/run1/asvtable_params_mfzr.tsv " \
              "--output asper2/run1/asvtable_params_taxa_mfzr.tsv --taxonomy vtam_db/taxonomy.tsv " \
              "--blastdbdir vtam_db/coi_blast_db --blastdbname coi_blast_db_20191211 -v --log asper2/vtam.log"

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=self.outdir_path)

    def test_04_zfzr_filter_optimized(self):

        cmd = "vtam filter --db asper2/db.sqlite --readinfo asper2/user_input/readinfo_zfzr.tsv " \
              "--readdir asper2/run1/sorted --params asper2/user_input/params_zfzr.yml " \
              "--asvtable asper2/run1/asvtable_params_zfzr.tsv -v --log asper2/vtam.log"

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=self.outdir_path)

    def test_05_zfzr_taxassign_optimized(self):

        cmd = "vtam taxassign --db asper2/db.sqlite --variants asper2/run1/asvtable_params_zfzr.tsv " \
              "--output asper2/run1/asvtable_params_taxa_zfzr.tsv --taxonomy vtam_db/taxonomy.tsv " \
              "--blastdbdir vtam_db/coi_blast_db --blastdbname coi_blast_db_20191211 -v --log asper2/vtam.log"

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=self.outdir_path)

    def test_07_pool(self):

        cmd = "vtam pool --db asper2/db.sqlite --runmarker asper2/user_input/pool_run_marker.tsv --output asper2/asvtable_pooled_mfzr_zfzr.tsv --log asper2/vtam.log -v"

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=self.outdir_path)

    def test_08_pool_taxa(self):

        cmd = "vtam taxassign --db asper2/db.sqlite --variants asper2/asvtable_pooled_mfzr_zfzr.tsv " \
              "--output asper2/asvtable_pooled_mfzr_zfzr_taxa.tsv --taxonomy vtam_db/taxonomy.tsv " \
              "--blastdbdir vtam_db/coi_blast_db --blastdbname coi_blast_db_20191211 --log asper2/vtam.log -v"

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=self.outdir_path)

    @classmethod
    def tearDownClass(cls):

        shutil.rmtree(cls.outdir_path, ignore_errors=True)
