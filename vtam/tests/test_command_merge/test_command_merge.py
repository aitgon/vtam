import filecmp
import os
import shutil
import unittest

from vtam import CommandMerge
from vtam.utils.PathManager import PathManager


class TestCommandMerge(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.outdir_path = os.path.join(PathManager.get_test_path(), 'outdir')

    def setUp(self):
        self.fastqinfo = os.path.join(os.path.dirname(__file__), "fastqinfo.tsv")
        self.fastqdir = os.path.dirname(__file__)

        self.fastainfo = os.path.join(self.outdir_path, "fastainfo.tsv")
        self.fastadir = os.path.join(self.outdir_path, "merged")

        self.fastainfo_bak = os.path.join(os.path.dirname(__file__), "fastainfo.tsv")
        self.fastadir_bak = self.outdir_path

    def test_01(self):

        CommandMerge.main(fastqinfo=self.fastqinfo, fastqdir=self.fastqdir, fastainfo=self.fastainfo,
                          fastadir=self.fastadir)
        self.assertTrue(filecmp.cmp(self.fastainfo, self.fastainfo_bak, shallow=True))
        self.assertTrue(filecmp.cmpfiles(self.fastadir, self.fastadir_bak, common=[
            'MFZR_14Ben01_1_fw_48.fasta', 'MFZR_14Ben01_2_fw_48.fasta'], shallow=True))

    def tearDown(self):
        shutil.rmtree(self.outdir_path, ignore_errors=True)
