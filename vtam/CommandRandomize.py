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
from vtam.utils.FileSampleInformation import FileSampleInformation


class CommandRandomize(object):
    """Class for the Randomize command"""
                
    @staticmethod
    def main(fastadir, random_seqdir, fastainfo, random_seqinfo, samplesize):

        if not os.path.isdir(fastadir) or not os.listdir(fastadir):
            Logger.instance().error(f"{fastadir} is empty or does not exists!")
            return

        fastainfo_df = FileSampleInformation(fastainfo).read_tsv_into_df()
        input_files = fastainfo_df.to_dict(orient='list')['mergedfasta']

        fastadir_path = os.path.abspath(fastadir)


        # check number is not > the sizes of the fasta files in fastadir 
        files_size = {}
        
        for input_file in input_files:
                
            file_path = os.path.join(fastadir_path, input_file)

            line_counter = LineCounter(file_path)
            file_size = line_counter.sequence_counter()

            files_size[input_file] = file_size

        smallest = min(files_size.values())

        if smallest < samplesize:
            Logger.instance().error(f"The smallest file in fastadir has {smallest} sequences.\nSamplesize cannot exceed this number of sequences")
            return


        ###################################################################
        #
        # Make the random files
        #
        ###################################################################
        
        # create output folder
        pathlib.Path(random_seqdir).mkdir(parents=True, exist_ok=True)

        output_files = []
        input_files_no_repeat = []

        # create the list to put in the ouput csv file and the list of file to sample (no duplicate)
        for input_file in input_files:
            base, ext = input_file.split(".", 1)
            output_files.append(base + "_sampled." + ext)
            if input_file not in input_files_no_repeat:
                input_files_no_repeat.append(input_file)


        for input_file in input_files_no_repeat:
            
            # get random indexes of lines in the file 
            lines = [] 
            num = 0
            for _ in range(samplesize):
                while num in lines:
                    num = randint(0, files_size[input_file])
                lines.append(num)

            lines = sorted(lines)
            
            # make output file path
            base, ext = input_file.split(".", 1)
            output_file = os.path.join(random_seqdir, base + "_sampled." + ext)

            # check extension
            if input_file.endswith(".gz"):
                _open = partial(gzip.open) 
            elif input_file.endswith(".bz2"):
                _open = partial(bz2.open)
            else:
                _open = open
            
            input_file = os.path.join(fastadir_path, input_file)
            
            with _open(input_file, 'rb') as f_in:
                with _open(output_file, 'ab') as f_out:
                    countline = -1
                    countSelect = 0
                    for line in f_in:
                        if line.startswith(b">"):
                            countline += 1

                            if countSelect + 1 < len(lines) and countline == lines[countSelect + 1]:
                                countSelect += 1
                        if countline == lines[countSelect]:
                            f_out.write(line)

                        if countline > lines[-1]:
                            break
            
        random_seqinfo_df = fastainfo_df.copy()
        random_seqinfo_df['mergedfasta'] = output_files
        random_seqinfo_df.to_csv(random_seqinfo, sep="\t", header=True, index=False)     
                       
                        




        