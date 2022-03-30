import os
import pathlib
import gzip 
import bz2
from functools import partial
from random import randint

from vtam.utils.VTAMexception import VTAMexception
from vtam.utils.Logger import Logger
from vtam.utils.FileCompression import FileCompression
from vtam.utils.LineCounter import LineCounter


class CommandRandomize(object):
    """Class for the Randomize command"""
                
    @staticmethod
    def main(fastadir, sampleddir, samplesize):

        input_files = os.listdir(fastadir)
        fastadir_path = os.path.abspath(fastadir)


        # check number is not > the sizes of the fasta files in fastadir 
        files_size = {}
        for input_file in input_files:

            file_path = os.path.join(fastadir_path, input_file)

            line_counter = LineCounter(file_path)
            file_size = line_counter.sequence_counter()

            if file_size < samplesize:
                Logger.instance().error(f"The largest file in fastadir has {file_size} sequences.\n samplesize cannot exceed this number of sequences")
                return
            files_size[input_file] = file_size

        print("files_size: \n", files_size)

        ###################################################################
        #
        # Make the random files
        #
        ###################################################################
        
        # creat output folder
        pathlib.Path(sampleddir).mkdir(parents=True, exist_ok=True)

        
        for input_file in input_files:
            
            # get random indexes of lines in the file 
            lines = [] #only unique odd numbers to be added
            num = 0
            for _ in range(samplesize):
                while num in lines:
                    num = randint(0, files_size[input_file])
                    if num % 2 != 0:
                        num += 1
                lines.append(num)
            lines = sorted(lines)
            print(f"lines for file: \n{lines}\n{input_file}")
            
            # make output file path
            base, ext = input_file.split(".", 1)
            output_file = os.path.join(sampleddir, base + "_sampled" + ext)
            
            # check extension
            if input_file.endswith(".gz"):      
                _open = partial(gzip.open) 
            elif input_file.endswith(".bz2"):
                _open = partial(bz2.open)
            else:
                _open = open
            
            input_file = os.path.join(fastadir_path, input_file)
            
            with open(input_file, 'rt') as f_in:
                with open(output_file, 'at') as f_out:
                    countline = -1
                    countSelect = 0
                    while countline <= lines[-1]:
                        line = f_in.readline()
                        if line.startswith(">"):
                            countline += 1

                            if countSelect + 1 < len(lines) and countline == lines[countSelect + 1]:
                                countSelect += 1
                        if countline == lines[countSelect]:
                            f_out.write(line)


                            
                       
                        




        