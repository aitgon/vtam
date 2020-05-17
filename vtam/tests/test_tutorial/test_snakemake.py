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
                 "Test requires online connection!")
class TestTutorialSnakemake(unittest.TestCase):

    """Will test main commands based on a complete test dataset"""

    @classmethod
    def setUpClass(cls):

        # vtam needs to be in the tsv_path
        subprocess.run([sys.executable, '-m', 'pip', 'install', '{}/.'.format(PathManager.get_package_path()),
                        '--upgrade'])

        cls.package_path = PathManager.get_package_path()
        cls.test_path = PathManager.get_test_path()
        cls.test_outdir_path = os.path.join(cls.test_path, 'outdir')
        shutil.rmtree(cls.test_outdir_path, ignore_errors=True)
        pathlib.Path(cls.test_outdir_path).mkdir(parents=True, exist_ok=True)

        # Set test paths
        cls.fastqinfo_path = os.path.join(cls.package_path, "doc/data/fastqinfo.tsv")
        cls.fastqdir_path = os.path.join(cls.test_outdir_path, "fastq")
        cls.snakefile_tuto_data = os.path.join(cls.package_path, "tools/snake.tuto.data.yml")
        cls.snakefile = os.path.join(cls.package_path, "doc/data/snakefile.yml")
        cls.snakeconfig_mfzr = os.path.join(cls.package_path, "doc/data/snakeconfig_mfzr.yml")

        ############################################################################################
        #
        # Download fastq test dataset
        #
        ############################################################################################

        fastq_tar_path = os.path.join(cls.test_outdir_path, "fastq.tar.gz")
        if not os.path.isfile(fastq_tar_path):
            urllib.request.urlretrieve(fastq_tar_gz_url, fastq_tar_path)
        tar = tarfile.open(fastq_tar_path, "r:gz")
        tar.extractall(path=cls.test_outdir_path)
        tar.close()

        cls.log_path = os.path.join(cls.test_outdir_path, "vtam.log")

        cls.asvtable_path = os.path.join(cls.test_outdir_path, "asvtable_default.tsv")

        cls.args = {}

        cls.args['package_path'] = cls.package_path
        cls.args['test_outdir_path'] = cls.test_outdir_path
        cls.args['fastqinfo'] = cls.fastqinfo_path
        cls.args['fastqdir'] = cls.fastqdir_path
        cls.args['snakefile'] = cls.snakefile
        cls.args['snakeconfig_mfzr'] = cls.snakeconfig_mfzr

        cls.args['snake_tuto_data'] = cls.snakefile_tuto_data

    def test_01(self):

        ############################################################################################
        #
        # Command Merge
        #
        ############################################################################################

        command = "snakemake --cores 1 -d {test_outdir_path} -s {snake_tuto_data} --config MARKER=mfzr " \
                  "PROJECT=asper1 PACKAGE_PATH={package_path} --until all_one_marker".format(**self.args)
        subprocess.run(shlex.split(command), check=True)

        command = "snakemake --directory {test_outdir_path} --resources db=1 " \
                    "--snakefile {snakefile} --cores 1 " \
                    "--configfile {snakeconfig_mfzr} --until asvtable_taxa".format(**self.args)
        subprocess.run(shlex.split(command), check=True)

        run1_mfzr_asvtable_default_taxa = os.path.join(self.test_outdir_path, "asper1/run1_mfzr/asvtable_default_taxa.tsv")
        run1_mfzr_asvtable_default_taxa_bak = os.path.join(self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_asvtable_default_taxa.tsv")
        self.assertTrue(filecmp.cmp(run1_mfzr_asvtable_default_taxa_bak, run1_mfzr_asvtable_default_taxa, shallow=True))

    @classmethod
    def tearDownClass(cls):

        subprocess.run([sys.executable, '-m', 'pip', 'uninstall', 'vtam', '-y'])
        shutil.rmtree(cls.test_outdir_path, ignore_errors=True)

