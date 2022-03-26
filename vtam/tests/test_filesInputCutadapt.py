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

    @classmethod
    def setUpClass(cls):
        cls.test_path = PathManager.get_test_path() 
        cls.fastainfo = os.path.join(cls.test_path, "test_files", "fastainfo.tsv")
        cls.mergedFasta1 = "MFZR_14Ben01_1_fw_48.fasta"
        cls.tags_file_path = os.path.join(cls.test_path, "test_files/FilesInputCutadapt")

    @staticmethod
    def compareFileContent(path1, path2):
        with open(path1, 'rt') as one:
            with open(path2, 'rt') as two:
                for line_one, line_two in zip(one.readlines(), two.readlines()):
                    if line_one != line_two:
                        return False
                return True

    # def setUp(self):
    #     self.fastainfo = os.path.join(self.test_path, "test_files", "fastainfo.tsv")


    # FilesInputCutadapt(self.fastainfo , self.mergedFasta1, no_reverse, tag_to_end, primer_to_end)
    def test_tags_file(self):
        #test with no option selected
        # os.path.join(self.outdir_path, "tags_file")
        filesCutadapt = FilesInputCutadapt(self.fastainfo , self.mergedFasta1, True, True, True)
        
        self.tags_file_path_result = os.path.join(self.test_path, filesCutadapt.tags_file())
        self.tags_file_path_ref = os.path.join(self.tags_file_path,"tagsFile1NoOption.fasta")
        #put the tags_file in a folder with a self path to be deleted
        self.assertTrue(self.compareFileContent(self.tags_file_path_result, self.tags_file_path_ref))

   
    def tearDown(self):
        if os.path.exists(self.tags_file_path_result):
            os.remove(self.tags_file_path_result)
        
