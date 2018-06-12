from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
import os
from wopmars.utils.Logger import Logger
import subprocess

class Merge(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.Merge"
    }

    __input_sample2fastq = "sample2fastq"
    __output_sample2fasta = "sample2fasta"

    def specify_input_file(self):
        return[Merge.__input_sample2fastq]

    def specify_output_file(self):
        return[Merge.__output_sample2fasta]

    def specify_params(self):
        return{
            "fastq_directory":"str",
            "fasta_dir":"str",
            "fastq_minovlen":"int",
            "fastq_maxmergelen":"int",
            "fastq_minmergelen":"int",
            "fastq_minlen":"int",
            "fastq_maxee":"int",
            "fastq_truncqual":"int",
            "fastq_maxns":"int",
            "threads":"int",
            "fastq_ascii":"int"
        }

    def run(self):
        # input files
        sample2fastq = self.input_file(Merge.__input_sample2fastq)
        # output files
        sample2fasta = self.output_file(Merge.__output_sample2fasta)

        try:
            # Opening the file to get all the lines stocked in the list csv_content
            with open(sample2fastq, 'r') as csv_file:
                next(csv_file)
                csv_content = csv_file.readlines()
            with open(sample2fasta, 'w') as csv_fout:
                for line in csv_content:
                    sample_info = line.strip().split("\t")
                    fin_fw = os.path.join(self.option("fastq_directory"), sample_info[8])
                    fin_rev = os.path.join(self.option("fastq_directory"), sample_info[9])
                    fout_name = sample_info[7] + "_" + sample_info[4] + "_" +sample_info[6] + ".fasta"
                    fout = os.path.join(self.option("fasta_dir"), fout_name)
                    line2write = line.strip() + '\t' + fout_name + "\n"
                    csv_fout.write(line2write)
                    if not os.path.isfile(fout):
                        if os.path.isfile(fin_fw) is False:
                            error_file = sample_info[8]
                            raise FileNotFoundError('One of the input file: ' + error_file + ' is missing.')

                        elif os.path.isfile(fin_rev) is False:
                            error_file = sample_info[9]
                            raise FileNotFoundError('One of the input file: ' + error_file + ' is missing.')
                        else:
                            command = "vsearch" + " -fastq_mergepairs " + fin_fw + " --reverse " + fin_rev
                            command = command + " --fastaout " + fout + ' --fastq_minovlen ' + str(self.option("fastq_minovlen"))
                            command = command + ' --fastq_maxmergelen ' + str(
                                self.option("fastq_maxmergelen")) + ' --fastq_minmergelen ' + str(self.option("fastq_minmergelen"))
                            command = command + " --fastq_minlen " + str(self.option("fastq_minlen")) + ' --fastq_maxee ' + str(
                                self.option("fastq_maxee"))
                            command = command + " --fastq_truncqual " + str(self.option("fastq_truncqual")) + " --fastq_maxns " + str(
                                self.option("fastq_maxns"))
                            command = command + " --threads " + str(self.option("threads")) + ' --fastq_ascii ' + str(
                                self.option("fastq_ascii"))
                            subprocess.call(command, shell=True)
                    else:
                        continue
        except IsADirectoryError:
            context = 'While trying to open csv file, the path is pointing to a directory and not a file'
            raise IsADirectoryError(context + ' Please check the path you enter for the csv option')
        except FileNotFoundError:
            context = \
                'While searching for the csv file or the fastq directory, any file or directory are found at the pointed directory'
            raise FileNotFoundError(context + ' Please check if the csv file or the fastq directory are present')
