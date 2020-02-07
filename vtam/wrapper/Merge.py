import sys

from wopmars.models.ToolWrapper import ToolWrapper
import os

import pathlib

from vtam.utils.VSearch import VSearch

from vtam import VTAMexception
from vtam.utils.Logger import Logger


class Merge(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.Merge"
    }

    __input_fastqinfo = "fastqinfo"
    __output_fastainfo = "fastainfo"

    def specify_input_file(self):
        return[Merge.__input_fastqinfo]

    def specify_output_file(self):
        return[Merge.__output_fastainfo]

    def specify_params(self):
        return{
            "fastq_dir": "str",
            "fasta_dir": "str",
            "fastq_minovlen": "int",
            "fastq_maxmergelen": "int",
            "fastq_minmergelen": "int",
            "fastq_minlen": "int",
            "fastq_maxee": "int",
            "fastq_truncqual": "int",
            "fastq_maxns": "int",
            "fastq_ascii": "int",
        }

    def run(self):
        session = self.session

        fastq_dir = str(os.getenv('VTAM_FASTQ_DIR'))
        fasta_dir = str(os.getenv('VTAM_FASTA_DIR'))

        ##########################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ##########################################################

        # input files
        fastqinfo = self.input_file(Merge.__input_fastqinfo)
        #
        # output files
        fastainfo = self.output_file(Merge.__output_fastainfo)
        #
        # Options
        fastq_minovlen = self.option("fastq_minovlen")
        fastq_maxmergelen = self.option("fastq_maxmergelen")
        fastq_minmergelen = self.option("fastq_minmergelen")
        fastq_minlen = self.option("fastq_minlen")
        fastq_maxee = self.option("fastq_maxee")
        fastq_truncqual = self.option("fastq_truncqual")
        fastq_maxns = self.option("fastq_maxns")
        fastq_ascii = self.option("fastq_ascii")
        #
        # Go
        # Opening the file to get all the lines stocked in the list csv_content
        fastq_and_fasta_list = [] # unique pairs of fastq files with the corresponding fasta_path file
        with open(fastqinfo, 'r') as csv_file:
            with open(fastainfo, 'w') as fastainfo_fout:
                next(csv_file) # skip header of fastqinfo
                fastainfo_header = "TagPair Forward	Primer Forward	TagPair Reverse	Primer Reverse	Marker name	 " \
                                   "Biosample	Replicate	Run	Fastq_fw	Fastq_rv	Fasta\n"
                fastainfo_fout.write(fastainfo_header) # write header of fastainfo
                for line in csv_file:
                    sample_info = line.strip().split("\t")
                    fastq_fw_abspath = os.path.join(fastq_dir, sample_info[8])
                    fastq_rv_abspath = os.path.join(fastq_dir, sample_info[9])
                    try:
                        pathlib.Path(fastq_fw_abspath).resolve(strict=True)
                    except FileNotFoundError:
                        Logger.instance().error(VTAMexception("VTAMexception: This FASTQ file was not found: {}.".format(fastq_fw_abspath)))
                        sys.exit(1)
                    try:
                        pathlib.Path(fastq_rv_abspath).resolve(strict=True)
                    except FileNotFoundError:
                        Logger.instance().error(VTAMexception("VTAMexception: This FASTQ file was not found: {}.".format(fastq_rv_abspath)))
                        sys.exit(1)
                    fasta_merged_basename = '.'.join(os.path.basename(fastq_fw_abspath).split('.')[0:-1]) + '_merged.fasta'
                    pathlib.Path(fasta_dir).mkdir(exist_ok=True)
                    fastainfo_line = line.strip() + '\t' + fasta_merged_basename + '\n'
                    fastainfo_fout.write(fastainfo_line) # write fastainfo line
                    fasta_abspath = os.path.join(fasta_dir, fasta_merged_basename)
                    if not (fastq_fw_abspath, fastq_rv_abspath, fasta_abspath) in fastq_and_fasta_list:
                        fastq_and_fasta_list.append((fastq_fw_abspath, fastq_rv_abspath, fasta_abspath))

        #########################################
        #
        # Loop of fastq pairs and run vsearch
        #
        #########################################
        for fastq_fw_abspath, fastq_rv_abspath, fasta_abspath in fastq_and_fasta_list:

            # Create object and run vsearch
            vsearch_parameters = {'--fastq_mergepairs': fastq_fw_abspath,
                                  '--reverse': fastq_rv_abspath,
                                  '--fastaout': fasta_abspath,
                                  '--fastq_minovlen': fastq_minovlen,
                                  '--fastq_maxmergelen': fastq_maxmergelen,
                                  '--fastq_minmergelen': fastq_minmergelen,
                                    "--fastq_minlen": fastq_minlen,
                                    "--fastq_maxee": fastq_maxee,
                                    "--fastq_truncqual": fastq_truncqual,
                                    "--fastq_maxns": fastq_maxns,
                                    "--fastq_ascii": fastq_ascii,
                                    "--threads": int(os.getenv('VTAM_THREADS')),
                                  }
            vsearch_cluster = VSearch(parameters=vsearch_parameters)
            vsearch_cluster.run()
