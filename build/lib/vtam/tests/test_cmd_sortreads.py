from unittest import TestCase
from vtam import CommandSortReads
from vtam.utils.PathManager import PathManager
import filecmp
import os
import shutil


class TestCommandSortReads(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.test_path = PathManager.get_test_path()
        cls.outdir_path = os.path.join(cls.test_path, 'outdir')

    def setUp(self):
        
        self.fastainfo = os.path.join(self.test_path, "test_files", "mergedinfo.tsv")
        self.fastadir = os.path.join(self.test_path, "test_files", "merged")
        self.sorted_dir = os.path.join(self.outdir_path, "sorted")
        self.sorted_dir_bak = os.path.join(self.test_path, "test_files", "sorted")

    def test_01(self):

        CommandSortReads.main(fastainfo=self.fastainfo, fastadir=self.fastadir, sorteddir=self.sorted_dir)

        self.assertTrue(filecmp.cmpfiles(self.sorted_dir, self.sorted_dir_bak, common=[
            'sortedinfo.tsv', 'MFZR_14Ben01_Tpos1_1_fw_48_000.fasta', 'MFZR_14Ben01_Tpos1_1_fw_48_001.fasta',
            'MFZR_14Ben01_Tpos1_1_fw_48_002.fasta', 'MFZR_14Ben01_Tpos1_1_fw_48_003.fasta'], shallow=True))

    def tearDown(self):
        shutil.rmtree(self.outdir_path, ignore_errors=True)
