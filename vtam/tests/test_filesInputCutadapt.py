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
        cls.tags_file_path = os.path.join(cls.test_path, "test_files", "FilesInputCutadapt")
        cls.fastainfo = os.path.join(cls.tags_file_path, "fastainfo.tsv")
        cls.mergedFasta1 = "14Ben01_1_fw_48.fasta"

    @staticmethod
    def compareFileContent(path1, path2):
        with open(path1, 'rt') as one:
            with open(path2, 'rt') as two:
                return one.read() == two.read()
                # for line_one, line_two in zip(one.readlines(), two.readlines()):
                #     if line_one != line_two:
                #         return False
                # return True

    def setUp(self):
        self.tags_file_path_result = ""


    def test_tags_file_noOption(self):
        # test with no option selected
        filesCutadapt = FilesInputCutadapt(self.fastainfo , self.mergedFasta1, True, True, True)
        self.tags_file_path_result = os.path.join(self.test_path, filesCutadapt.tags_file())
        self.tags_file_path_ref = os.path.join(self.tags_file_path,"tagsFile1NoOption.fasta")
        
        self.assertTrue(self.compareFileContent(self.tags_file_path_result, self.tags_file_path_ref))

    def test_tags_file_noReverse(self):
        # test with no_reverse selected
        filesCutadapt = FilesInputCutadapt(self.fastainfo , self.mergedFasta1, False, True, True)
        
        self.tags_file_path_result = os.path.join(self.test_path, filesCutadapt.tags_file())
        self.tags_file_path_ref = os.path.join(self.tags_file_path,"tagsFile1NoReverse.fasta")
        
        self.assertTrue(self.compareFileContent(self.tags_file_path_result, self.tags_file_path_ref))

    def test_tags_file_tagToEnd(self):
        # test with no_reverse selected
        filesCutadapt = FilesInputCutadapt(self.fastainfo , self.mergedFasta1, True, False, True)
        self.tags_file_path_result = os.path.join(self.test_path, filesCutadapt.tags_file())
        self.tags_file_path_ref = os.path.join(self.tags_file_path,"tagsFile1TagToEnd.fasta")
        
        self.assertTrue(self.compareFileContent(self.tags_file_path_result, self.tags_file_path_ref))

    def test_primers(self):
        
        primers_test_no_duplicate = [('mfzr', 'TCCACTAATCACAARGATATTGGTAC', 'WACTAATCAATTWCCAAATCCTCC', '26', '24'), ('zfzr', 'AGATATTGGAACWTTATATTTTATTTTTGG', 'WACTAATCAATTWCCAAATCCTCC', '30', '24')]
        primers_test_duplicate = [primers_test_no_duplicate[0]] + primers_test_no_duplicate

        filesCutadapt = FilesInputCutadapt(self.fastainfo , self.mergedFasta1, False, False, False)
        primers = filesCutadapt.primers()

        self.assertTrue(primers == primers_test_no_duplicate)
        self.assertFalse(primers == primers_test_duplicate)

    def test_get_df_info(self):

        dict_test = {
            'run': ['run1', 'run1', 'run1'],
            'marker': ['mfzr', 'mfzr', 'zfzr'],
            'sample': ['14ben01', '14ben01', '14ben01'],
            'replicate': [1, 1, 1],
            'tagfwd': ['gtcgatcatgtca', 'gtcgatcatgtca', 'gtcgatcatgtca'],
            'primerfwd': ['TCCACTAATCACAARGATATTGGTAC', 'TCCACTAATCACAARGATATTGGTAC', 'AGATATTGGAACWTTATATTTTATTTTTGG'],
            'tagrev': ['acatcgacgtacg', 'acatcgacgtacg', 'acatcgacgtacg'],
            'primerrev': ['WACTAATCAATTWCCAAATCCTCC', 'WACTAATCAATTWCCAAATCCTCC', 'WACTAATCAATTWCCAAATCCTCC'],
            'mergedfasta': ['14Ben01_1_fw_48.fasta', '14Ben01_1_fw_48.fasta', '14Ben01_1_fw_48.fasta']
        }

        filesCutadapt = FilesInputCutadapt(self.fastainfo , self.mergedFasta1, False, False, False)

        self.assertTrue(dict_test == filesCutadapt.get_df_info())
    
    def test_get_sample_names(self):
        names_test = ['14ben01']
        names_test_reversed = ['14ben01', '14ben01_reversed']

        filesCutadapt = FilesInputCutadapt(self.fastainfo , self.mergedFasta1, False, False, False)
        names = filesCutadapt.get_sample_names()
        self.assertTrue(names == names_test)

        filesCutadapt = FilesInputCutadapt(self.fastainfo , self.mergedFasta1, True, False, False)
        names = filesCutadapt.get_sample_names()
        self.assertTrue(names == names_test_reversed)


    def tearDown(self):
        if os.path.exists(self.tags_file_path_result):
            os.remove(self.tags_file_path_result)
        
