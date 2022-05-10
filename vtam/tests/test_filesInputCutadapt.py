import unittest 
import os

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
        cls.fastainfoNoDuplicates = os.path.join(cls.tags_file_path, "fastainfoNoDuplicates.tsv")
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
        filesCutadapt = FilesInputCutadapt(self.fastainfoNoDuplicates , self.mergedFasta1, True, True)
        self.tags_file_path_result = os.path.join(self.test_path, filesCutadapt.tags_file())
        self.tags_file_path_ref = os.path.join(self.tags_file_path,"tagsFile1NoOption.fasta")
        
        self.assertTrue(self.compareFileContent(self.tags_file_path_result, self.tags_file_path_ref))

    def test_tags_file_duplicates(self):
        # test with no option selected            
        filesCutadapt = FilesInputCutadapt(self.fastainfo , self.mergedFasta1, True, True)
        
        exception = "Exception: ('run[run1]marker[zfzr]sample[14ben01]replicate[1]fwd[gtcgatcatgtca]rev[acatcgacgtacg]', 'run1', 'zfzr', '14ben01', '1', 'gtcgatcatgtca', 'acatcgacgtacg') and ('run[run1]marker[mfzr]sample[14ben01]replicate[1]fwd[gtcgatcatgtca]rev[acatcgacgtacg]', 'run1', 'mfzr', '14ben01', '1', 'gtcgatcatgtca', 'acatcgacgtacg') lines have different run/marker/sample/replicate combinations but same tag_combination')"
        self.assertRaises(Exception, filesCutadapt.tags_file)

    def test_tags_file_noReverse(self):
        # test with no_reverse selected
        filesCutadapt = FilesInputCutadapt(self.fastainfoNoDuplicates , self.mergedFasta1, False, True)
        
        self.tags_file_path_result = os.path.join(self.test_path, filesCutadapt.tags_file())
        self.tags_file_path_ref = os.path.join(self.tags_file_path,"tagsFile1NoReverse.fasta")

        self.assertTrue(self.compareFileContent(self.tags_file_path_result, self.tags_file_path_ref))

    def test_tags_file_tagToEnd(self):
        # test with no_reverse selected
        filesCutadapt = FilesInputCutadapt(self.fastainfoNoDuplicates , self.mergedFasta1, True, False)
        self.tags_file_path_result = os.path.join(self.test_path, filesCutadapt.tags_file())
        self.tags_file_path_ref = os.path.join(self.tags_file_path,"tagsFile1TagToEnd.fasta")
        
        self.assertTrue(self.compareFileContent(self.tags_file_path_result, self.tags_file_path_ref))

    def test_primers(self):
        
        primers_test_no_duplicate = [('mfzr', 'TCCACTAATCACAARGATATTGGTAC', 'WACTAATCAATTWCCAAATCCTCC', 26, 24), ('zfzr', 'AGATATTGGAACWTTATATTTTATTTTTGG', 'WACTAATCAATTWCCAAATCCTCC', 30, 24)]
        primers_test_duplicate = [primers_test_no_duplicate[0]] + primers_test_no_duplicate

        filesCutadapt = FilesInputCutadapt(self.fastainfo , self.mergedFasta1, False, False)
        primers = filesCutadapt.primers()

        self.assertTrue(primers == primers_test_no_duplicate)
        self.assertFalse(primers == primers_test_duplicate)

    def test_get_df_info(self):

        dict_test = {
            'run': ['run1', 'run1'],
            'marker': ['mfzr', 'zfzr'],
            'sample': ['14ben01', '14ben01'],
            'replicate': [1, 1],
            'tagfwd': ['gtcgatcatgtca', 'gtcgatcatgtca'],
            'primerfwd': ['TCCACTAATCACAARGATATTGGTAC', 'AGATATTGGAACWTTATATTTTATTTTTGG'],
            'tagrev': ['acatcgacgtacg', 'acatcgacgtacg'],
            'primerrev': ['WACTAATCAATTWCCAAATCCTCC', 'WACTAATCAATTWCCAAATCCTCC'],
            'mergedfasta': ['14Ben01_1_fw_48.fasta','14Ben01_1_fw_48.fasta']
        }

        

        filesCutadapt = FilesInputCutadapt(self.fastainfo , '14Ben01_1_fw_48.fasta', False, False)
        self.assertTrue(dict_test == filesCutadapt.get_df_info())
    
    def test_get_sample_names(self):
        names_test = [('run[run1]marker[mfzr]sample[14ben01]replicate[1]fwd[gtcgatcatgtca]rev[acatcgacgtacg]', 'run1', 'mfzr', '14ben01', '1', 'gtcgatcatgtca', 'acatcgacgtacg')]
        names_test_reversed = [('run[run1]marker[mfzr]sample[14ben01]replicate[1]fwd[gtcgatcatgtca]rev[acatcgacgtacg]', 'run1', 'mfzr', '14ben01', '1', 'gtcgatcatgtca', 'acatcgacgtacg'), 
        ('run[run1]marker[mfzr]sample[14ben01]replicate[1]fwd[gtcgatcatgtca]rev[acatcgacgtacg]_reversed', 'run1', 'mfzr', '14ben01', '1', 'gtcgatcatgtca', 'acatcgacgtacg')]

        filesCutadapt = FilesInputCutadapt(self.fastainfoNoDuplicates , self.mergedFasta1, False, False)
        names = filesCutadapt.get_sample_names()
        self.assertTrue(names == names_test)

        filesCutadapt = FilesInputCutadapt(self.fastainfoNoDuplicates , self.mergedFasta1, True, False)
        names = filesCutadapt.get_sample_names()
        self.assertTrue(names == names_test_reversed)


    def tearDown(self):
        if os.path.exists(self.tags_file_path_result):
            os.remove(self.tags_file_path_result)
        
