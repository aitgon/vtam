
import shutil
from unittest import TestCase

from vtam import CommandSortReads
from vtam.utils.PathManager import PathManager

import os
import filecmp


class TestCommandSortReads(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.outdir_path = os.path.join(PathManager.get_test_path(), 'outdir')

    def setUp(self):
        self.fastainfo = os.path.join(os.path.dirname(__file__), "fastainfo.tsv")
        self.fastadir = os.path.join(os.path.dirname(__file__), "merged")

        self.sorted_dir = os.path.join(self.outdir_path, "sorted")

        self.sorted_dir_bak = os.path.join(self.outdir_path, 'sorted')

    def test_01(self):

        CommandSortReads.main(fastainfo=self.fastainfo, fastadir=self.fastadir, outdir=self.sorted_dir)

        self.assertTrue(filecmp.cmpfiles(self.sorted_dir, self.sorted_dir_bak, common=[
            'readinfo.tsv', 'MFZR_14Ben01_Tpos1_1_fw_48_000.fasta', 'MFZR_14Ben01_Tpos1_1_fw_48_001.fasta',
            'MFZR_14Ben01_Tpos1_1_fw_48_002.fasta', 'MFZR_14Ben01_Tpos1_1_fw_48_003.fasta'], shallow=True))

    def tearDown(self):
        shutil.rmtree(self.outdir_path, ignore_errors=True)
