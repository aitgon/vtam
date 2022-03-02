import unittest 
import os
from shutil import copyfile
import filecmp

from vtam.utils.PathManager import PathManager
from vtam.utils.FileDecompression import FileDecompression
# import imp    
# FileDecompression = imp.load_source('FileDecompression', '../../vtam/utils/py')

class TestFileDecompression(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.test_path = PathManager.get_test_path() # return the path vtam.test_path__path__[0]/tests
        cls.outdir_path = os.path.join(cls.test_path, 'outdir')

    def setUp(self):

        # PATHS to dir
        self.fastq_dir = os.path.join(self.test_path, "test_files", "fastq")
        self.fastq_gz_dir = os.path.join(self.test_path, "test_files", "fastq_gz")
        # PATHS to files
        self.fastq_gz_file = os.path.join(self.fastq_gz_dir, "MFZR_14Ben01_1_fw_48.fastq.gz")
        self.fastq_pigz_file = os.path.join(self.fastq_gz_dir, "MFZR_14Ben01_1_fw_48_pigz.fastq.gz")
        self.fastq_file = os.path.join(self.fastq_dir, "MFZR_14Ben01_1_fw_48.fastq")

        #TODO: add tests for .bz2 files

        # Make a copy of the file to be deleted
        self.fastq_gz_file_copy = os.path.join(self.fastq_gz_dir, "MFZR_14Ben01_1_fw_48_copy.fastq.gz")
        copyfile(self.fastq_gz_file, self.fastq_gz_file_copy)


    def test_gzip_decompression(self):

        decompression_instance = FileDecompression(self.fastq_gz_file_copy)
        self.decompressed = decompression_instance.gzip_decompression()

        ## test if the returned file name has the correct extension
        self.assertTrue(self.decompressed.endswith('.fastq'))
        ## test filename of the result when given a file which is already .gz format
        self.assertFalse(self.decompressed.endswith('.fastq.gz'))
        ## test if the file content of the returned file is compressed
        self.assertTrue(filecmp.cmp(self.decompressed, self.fastq_file, shallow=True))
        
        decompression_instance.delete_file()
        #test if original file gets deleted
        self.assertFalse(os.path.exists(self.fastq_gz_file_copy))

        #test if wrong path is given
        decompression_instance_wrong = FileDecompression("wrong/file/path")
        compressed_wrong = decompression_instance_wrong.gzip_decompression()
        self.assertFalse(compressed_wrong)
        self.assertFalse(os.path.exists(self.fastq_gz_file_copy))


    def test_pigz_decompression(self):

        decompression_instance = FileDecompression(self.fastq_gz_file_copy)
        self.decompressed = decompression_instance.pigz_decompression()

        ## test if the returned file name has the correct extension
        self.assertFalse(self.decompressed.endswith('.gz'))        ## test filename of the result when given a file which is already .gz format
        self.assertTrue(self.decompressed.endswith('.fastq'))
        ## test if the file content of the returned file is compressed
        self.assertTrue(filecmp.cmp(self.decompressed, self.fastq_file, shallow=True))
        
        decompression_instance.delete_file()
        #test if original file gets deleted
        self.assertFalse(os.path.exists(self.fastq_gz_file_copy))

        #test if wrong path is given
        decompression_instance_wrong = FileDecompression("wrong/file/path")
        compressed_wrong = decompression_instance_wrong.pigz_decompression()
        self.assertFalse(compressed_wrong)


    def tearDown(self):
        if os.path.exists(self.fastq_gz_file_copy):
            os.remove(self.fastq_gz_file_copy)
        if os.path.exists(self.decompressed):
            os.remove(self.decompressed)

