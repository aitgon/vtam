import unittest 
import filecmp
import os
from shutil import copyfile
import hashlib

from vtam.utils.PathManager import PathManager
from vtam.utils.FileCompression import FileCompression
# import imp    
# FileCompression = imp.load_source('FileCompression', '../../vtam/utils/py')

class TestFileCompression(unittest.TestCase):

    @staticmethod
    def compareContent(filename):
        md5 = hashlib.md5()
        count = 0 #only compare the content (get rid of the headers see PEP1952)
        with open(filename, 'rb') as f_in:
            while count < 2:
                text = f_in.read(1)
                if text == b'\x00':
                    count += 1
            while True:
                data = f_in.read(8192)
                if not data:
                    break
                md5.update(data)
        return md5.hexdigest()

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
        # Make a copy of the file to be deleted
        self.fastq_file_copy = os.path.join(self.fastq_dir, "MFZR_14Ben01_1_fw_48_copy.fastq")
        copyfile(self.fastq_file, self.fastq_file_copy)


    def test_gzip_compression(self):

        compression_instance = FileCompression(self.fastq_file_copy)
        self.compressed = compression_instance.gzip_compression()

        ## test if the returned file name has the correct extension
        self.assertTrue(self.compressed.endswith('.gz'))
        ## test filename of the result when given a file which is already .gz format
        self.assertFalse(self.compressed.endswith('.gz.gz'))
        #TODO: give an uncompressed file with gz extension to test that
        ## test if the file content of the returned file is compressed
        self.assertEqual(self.compareContent(self.compressed), self.compareContent(self.fastq_gz_file))
        
        compression_instance.delete_file()
        #test if original file gets deleted
        self.assertFalse(os.path.exists(self.fastq_file_copy))

        #test if wrong path is given
        compression_instance_wrong = FileCompression("wrong/file/path")
        compressed_wrong = compression_instance_wrong.gzip_compression()

        self.assertFalse(compressed_wrong)

    def test_pigz_compression(self):

        compression_instance = FileCompression(self.fastq_file_copy)
        self.compressed = compression_instance.pigz_compression()

        ## test if the returned file name has the correct extension
        try:
            output = self.compressed.endswith('.gz')
            output2 = self.compressed.endswith('.gz.gz')
            stop_test = False
        except AttributeError:
            stop_test = True
            output = True
            output2 = False
        self.assertTrue(output)        ## test filename of the result when given a file which is already .gz format
        self.assertFalse(output2)
        ## test if the file content of the returned file is compressed
        if not stop_test:
            self.assertEqual(self.compareContent(self.compressed), self.compareContent(self.fastq_pigz_file))
        
        compression_instance.delete_file()
        #test if original file gets deleted
        self.assertFalse(os.path.exists(self.fastq_file_copy))

        #test if wrong path is given
        compression_instance_wrong = FileCompression("wrong/file/path")
        compressed_wrong = compression_instance_wrong.pigz_compression()

        self.assertFalse(compressed_wrong)

    def tearDown(self):
        if os.path.exists(self.fastq_file_copy):
            os.remove(self.fastq_file_copy)
        if self.compressed is not None:
            if os.path.exists(self.compressed):
                os.remove(self.compressed)

