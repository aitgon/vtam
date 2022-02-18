import filecmp
import os
import shutil
import unittest

from vtam import CommandMerge
from vtam.utils.PathManager import PathManager


class TestCommandMerge(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.test_path = PathManager.get_test_path() # return the path vtam.test_path__path__[0]/tests
        cls.outdir_path = os.path.join(cls.test_path, 'outdir_gz')

    def setUp(self):

        self.fastqinfo = os.path.join(self.test_path, "test_files", "fastqinfo_gz.tsv")
        self.fastqdir = os.path.join(self.test_path, "test_files", "fastq_gz")

        self.fastainfo = os.path.join(self.outdir_path, "fastainfo_gz.tsv")
        self.fastadir = os.path.join(self.outdir_path, "merged_gz")

        self.fastainfo_bak = os.path.join(self.test_path, "test_files", "fastainfo_gz.tsv")
        self.fastadir_bak = self.outdir_path

    def test_01(self):

        CommandMerge.main(fastqinfo=self.fastqinfo, fastqdir=self.fastqdir, fastainfo=self.fastainfo,
                          fastadir=self.fastadir)
        import pdb; pdb.set_trace()
        self.assertTrue(filecmp.cmp(self.fastainfo, self.fastainfo_bak, shallow=True))
        self.assertTrue(filecmp.cmpfiles(self.fastadir, self.fastadir_bak, common=[
            'MFZR_14Ben01_1_fw_48.fasta.gz', 'MFZR_14Ben01_2_fw_48.fasta'], shallow=True))

    def tearDown(self):
        shutil.rmtree(self.outdir_path, ignore_errors=True)
