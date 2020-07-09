# -*- coding: utf-8 -*-
import shlex
import sys
import filecmp
import os
import pathlib
import shutil
import subprocess
import tarfile
import urllib
import unittest

from vtam.utils.constants import fastq_tar_gz_url
from vtam.utils.PathManager import PathManager
from urllib import request


@unittest.skipIf(request.urlopen(fastq_tar_gz_url).getcode() != 200,
                 "This test requires an internet connection!")
class TestTutorialCommands(unittest.TestCase):

    """Will test main commands based on a complete test dataset"""

    @classmethod
    def setUpClass(cls):

        cls.package_path = os.path.join(PathManager.get_package_path())

        # vtam needs to be in the tsv_path
        subprocess.run([sys.executable, '-m', 'pip', 'install', os.path.join('{}'.format(cls.package_path), '.'),
                        '--upgrade'])

        cls.test_path = os.path.join(PathManager.get_test_path())
        cls.outdir_path = os.path.join(cls.test_path, 'outdir')
        shutil.rmtree(cls.outdir_path, ignore_errors=True)  # during development of the test, this prevents errors
        pathlib.Path(cls.outdir_path).mkdir(parents=True, exist_ok=True)

        ################################################################################################################
        #
        # Download fastq test dataset
        #
        ################################################################################################################

        fastq_tar_path = os.path.join(cls.outdir_path, "fastq.tar.gz")
        if not os.path.isfile(fastq_tar_path):
            urllib.request.urlretrieve(fastq_tar_gz_url, fastq_tar_path)
        tar = tarfile.open(fastq_tar_path, "r:gz")
        tar.extractall(path=cls.outdir_path)
        tar.close()

        # Set test paths
        cls.fastqinfo_path = os.path.join(PathManager.get_package_path(), "doc", "data", "fastqinfo.tsv")
        cls.fastqdir_path = os.path.join(cls.outdir_path, "fastq")
        cls.fastainfo_path = os.path.join(cls.outdir_path, "fastainfo.tsv")
        cls.fastadir_path = os.path.join(cls.outdir_path, "merged")

        cls.sorted_dir_path = os.path.join(cls.outdir_path, "sorted")
        cls.sortedreadinfo_path = os.path.join(cls.sorted_dir_path, "readinfo.tsv")

        cls.log_path = os.path.join(cls.outdir_path, "vtam.log")

        cls.asvtable_path = os.path.join(cls.outdir_path, "asvtable_default.tsv")

        cls.args = {}
        cls.args['fastqinfo'] = cls.fastqinfo_path
        cls.args['fastqdir'] = cls.fastqdir_path
        cls.args['fastainfo'] = cls.fastainfo_path
        cls.args['fastadir'] = cls.fastadir_path
        cls.args['sorted'] = cls.sorted_dir_path
        cls.args['db'] = os.path.join(cls.outdir_path, "db.sqlite")
        cls.args['readinfo'] = cls.sortedreadinfo_path
        cls.args['readdir'] = cls.sorted_dir_path
        cls.args['asvtable'] = cls.asvtable_path
        cls.args['log'] = cls.log_path

    def test_step01_merge(self):

        ################################################################################################################
        #
        # Command Merge
        #
        ################################################################################################################

        cmd = "vtam merge --fastqinfo {fastqinfo} --fastqdir {fastqdir} --fastainfo {fastainfo} --fastadir {fastadir} " \
              "-v --log {log}".format(**self.args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args)

        self.fastainfo_path_bak = os.path.join(self.test_path, "test_files_dryad.f40v5_small", "run1_mfzr_zfzr", "fastainfo.tsv")
        self.fastadir_path_bak = os.path.join(os.path.dirname(__file__), "merge")

        self.assertTrue(filecmp.cmp(self.fastainfo_path, self.fastainfo_path_bak, shallow=True))
        self.assertTrue(os.path.getsize(os.path.join(self.fastadir_path, 'mfzr_1_fw.fasta')) >= 11608260)
        self.assertTrue(os.path.getsize(os.path.join(self.fastadir_path, 'mfzr_1_fw.fasta')) <= 11795030)
        self.assertTrue(os.path.getsize(os.path.join(self.fastadir_path, 'zfzr_3_fw.fasta')) >= 11658700)
        self.assertTrue(os.path.getsize(os.path.join(self.fastadir_path, 'zfzr_3_fw.fasta')) <= 11838710)

    def test_step02_sortreads(self):

        ################################################################################################################
        #
        # Command SortReads
        #
        ################################################################################################################

        cmd = "vtam sortreads --fastainfo {fastainfo} --fastadir {fastadir} --outdir {sorted} " \
              "-v --log {log}".format(**self.args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args)

        self.sortedreadinfo_path_bak = os.path.join(self.test_path, "test_files_dryad.f40v5_small", "run1_mfzr_zfzr", "sortedreadinfo.tsv")
        self.assertTrue(filecmp.cmp(self.sortedreadinfo_path, self.sortedreadinfo_path_bak, shallow=True))
        self.assertTrue(os.path.getsize(os.path.join(self.sorted_dir_path, 'mfzr_1_fw_000.fasta')) >= 5131890)  # 5131896
        self.assertTrue(os.path.getsize(os.path.join(self.sorted_dir_path, 'mfzr_1_fw_000.fasta')) <= 5131900)
        self.assertTrue(os.path.getsize(os.path.join(self.sorted_dir_path, 'zfzr_3_fw_023.fasta')) >= 909500)  # 909507
        self.assertTrue(os.path.getsize(os.path.join(self.sorted_dir_path, 'zfzr_3_fw_023.fasta')) <= 909510)

    @classmethod
    def tearDownClass(cls):

        shutil.rmtree(cls.outdir_path, ignore_errors=True)

