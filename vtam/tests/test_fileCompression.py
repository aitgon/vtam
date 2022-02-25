import unittest 
import filecmp
import os
from shutil import copyfile
import hashlib

from vtam.utils.PathManager import PathManager
import imp    
FileCompression = imp.load_source('FileCompression', '../../vtam/utils/FileCompression.py')

class TestFileCompression(unittest.TestCase):

    @staticmethod
    def compareContent(filename):
        md5 = hashlib.md5()
        with open(filename, 'rb') as f_in:
            while True:
                text = f_in.read(1)
                if text == b'\x00':
                    break
            while True:
                data = f_in.read(8192)
                if not data:
                    break
                md5.update(data)
        return md5.hexdigest()

    @classmethod
    def setUpClass(cls):

        cls.test_path = os.getcwd()#PathManager.get_test_path() # return the path vtam.test_path__path__[0]/tests
        cls.outdir_path = os.path.join(cls.test_path, 'outdir')

    def setUp(self):

        # PATHS to dir
        self.fastq_dir = os.path.join(self.test_path, "test_files", "fastq")
        self.fastq_gz_dir = os.path.join(self.test_path, "test_files", "fastq_gz")
        # PATHS to files
        self.fastq_gz_file = os.path.join(self.fastq_gz_dir, "MFZR_14Ben01_1_fw_48.fastq.gz")
        self.fastq_file = os.path.join(self.fastq_dir, "MFZR_14Ben01_1_fw_48.fastq")
        # Make a copy of the file to be deleted
        self.fastq_gz_file_copy = os.path.join(self.fastq_gz_dir, "MFZR_14Ben01_1_fw_48_copy.fastq.gz")
        self.fastq_file_copy = os.path.join(self.fastq_dir, "MFZR_14Ben01_1_fw_48_copy.fastq")
        copyfile(self.fastq_file, self.fastq_file_copy)
        copyfile(self.fastq_gz_file, self.fastq_gz_file_copy) # to use if the filename is considered when filecmp


    def test_FileCompression(self):

        compression_instance = FileCompression.FileCompression(self.fastq_file_copy)
        self.compressed = compression_instance.gzip_compression()

        #import pdb; pdb.set_trace()

        ## test if the returned file name has the correct extension
        self.assertEqual(self.compressed, str(self.fastq_file_copy) + ".gz")
        ## test filename of the result when given a file which is already .gz format
        self.assertFalse(self.compressed.endswith('.gz.gz'))
        ## test if the file content of the returned file is compressed
        print(self.compressed)
        print(self.fastq_gz_file_copy)
        self.assertEqual(self.compareContent(self.compressed), self.compareContent(self.fastq_gz_file_copy))
        

        compression_instance.delete_file()

        #test if original file gets deleted
        self.assertFalse(os.path.exists(self.fastq_file_copy))

        #test if wrong path is given
        compression_instance_wrong = FileCompression.FileCompression("wrong/file/path")
        compressed_wrong = compression_instance_wrong.gzip_compression()

        self.assertFalse(compression_instance_wrong.file_path)

    def tearDown(self):
        if os.path.exists(self.fastq_file_copy):
            os.remove(self.fastq_file_copy)
        if os.path.exists(self.compressed):
            os.remove(self.compressed)

if __name__=="__main__":
    unittest.main()