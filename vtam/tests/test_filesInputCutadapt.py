import unittest 
import filecmp
import os
from shutil import copyfile
import hashlib

from vtam.utils.PathManager import PathManager
from vtam.utils.FilesInputCutadapt import FilesInputCutadapt
# import imp    
# FileCompression = imp.load_source('FileCompression', '../../vtam/utils/py')

class TestFilesInputCutadapt(unittest.TestCase):

    def setUpClass(cls):
        cls.test_path = PathManager.get_test_path()
        cls.outdir_path = os.path.join(cls.test_path, 'outdir')
        cls.fastainfo = os.path.join(cls.test_path, "test_files", "fastainfo.tsv")
        cls.mergedFasta1 = "MFZR_14Ben01_1_fw_48.fasta"

 


    # def setUp(self):
    #     self.fastainfo = os.path.join(self.test_path, "test_files", "fastainfo.tsv")


    # FilesInputCutadapt(self.fastainfo , self.mergedFasta1, no_reverse, tag_to_end, primer_to_end)
    def test_tags_file(self):
        #test with no option selected
        filesCutadapt = FilesInputCutadapt(self.fastainfo , self.mergedFasta1, True, True, True)
        tags_file = filesCutadapt.tags_file()

        #put the tags_file in a folder with a self path to be deleted
        self.assertTrue(filecmp.cmp(tags_file, self.tags_file, shallow=True))






        
        

   
    def tearDown(self):
        if os.path.exists(self.fastq_file_copy):
            os.remove(self.fastq_file_copy)
        if self.compressed is not None:
            if os.path.exists(self.compressed):
                os.remove(self.compressed)

