import gzip
import bz2
import shutil
import os 
import shlex
import sys
import subprocess


class FileCompression(object):
    """ compress file to desired format """

    def __init__(self, file_path):
        self.file_path = file_path
        if not os.path.exists(self.file_path):
            self.file_path = False

        self.compressed = None

    # def find_files(path_to_folder):
    #     ''' iterates over a folder content and return a list of its content '''

    #     list_of_file = listdir(path_to_folder)

    #     return list_of_file

    #TODO add checks if file exists, read file extension (may be already compressed to another format) 

    def gzip_compression(self):
        ''' take a file and compress it to gz format using the gzip package'''

        if not self.file_path:
            return self.file_path

        if self.file_path.endswith('.gz'):
            path_to_compressed_file = self.file_path
        else:
            path_to_compressed_file = self.file_path + '.gz'

        with open(self.file_path, 'rb') as f_in:
            with gzip.open(path_to_compressed_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        if os.path.exists(path_to_compressed_file):
            self.compressed = path_to_compressed_file
            return path_to_compressed_file
        return self.file_path
    
    def pigz_compression(self):

        if not self.file_path:
            return self.file_path

        if self.file_path.endswith(".gz"):
            self.file_path = self.file_path[:-3]

        cmd = "pigz " + self.file_path

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        try :
            subprocess.run(args=args, check=True)
            self.compressed = self.file_path + ".gz"
            return self.compressed
        except:
            self.gzip_compression()
    
    def bz2_compression(self):
        ''' take a file and compress it to .bz2 format using the bz2 package'''

        if not self.file_path:
            return self.file_path

        if self.file_path.endswith('.bz2'):
            path_to_compressed_file = self.file_path
        else:
            path_to_compressed_file = self.file_path + '.bz2'

        with open(self.file_path, 'rb') as f_in:
            with bz2.open(path_to_compressed_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        if os.path.exists(path_to_compressed_file):
            self.compressed = path_to_compressed_file
            return path_to_compressed_file
        return self.file_path
    
    def delete_file(self):
        if self.file_path and self.compressed != self.file_path:
            os.remove(self.file_path)
        return self.file_path
         
