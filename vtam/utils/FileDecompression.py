import gzip
import shutil
import os 

class FileDecompression(object):
    """ Decompress file to from .gz format """

    def __init__(self, file_path):
        self.file_path = file_path
        if not os.path.exists(self.file_path):
            self.file_path = False

        self.decompressed = None

    #TODO add checks if file exists, read file extension (may be already compressed to another format) 

    def gzip_decompression(self):
        """ take a file and decompress it to .gz format using the gzip package """

        if not self.file_path:
            return self.file_path

        if self.file_path.endswith('.gz'):
            path_to_decompressed_file = self.file_path[:-3]
        else:
            return(self.path)

        with gzip.open(self.file_path, 'rb') as f_in:
            with open(path_to_decompressed_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        if os.path.exists(path_to_decompressed_file):
            self.decompressed = path_to_decompressed_file
            return self.decompressed
        return self.file_path

    def delete_file(self):
        if self.file_path and self.decompressed != self.file_path:
            os.remove(self.file_path)
        return self.file_path
         
