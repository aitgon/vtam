import sys

from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
import os

import subprocess
import pathlib

from wopmetabarcoding import VTAMexception
from wopmetabarcoding.utils.Logger import Logger
from wopmetabarcoding.utils.OptionManager import OptionManager
from wopmetabarcoding.utils.PathManager import PathManager


class Merge(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.Merge"
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
            "threads": "int",
            "fastq_ascii": "int",
            "log_verbosity": "int",
            "log_file": "str"
        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        OptionManager.instance()['log_verbosity'] = int(self.option("log_verbosity"))
        OptionManager.instance()['log_file'] = str(self.option("log_file"))

        ##########################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ##########################################################
        #
        # input files
        fastqinfo = self.input_file(Merge.__input_fastqinfo)
        #
        # output files
        fastainfo = self.output_file(Merge.__output_fastainfo)
        #
        # Options
        fastq_dir = self.option("fastq_dir")
        fasta_dir = self.option("fasta_dir")
        fastq_minovlen = self.option("fastq_minovlen")
        fastq_maxmergelen = self.option("fastq_maxmergelen")
        fastq_minmergelen = self.option("fastq_minmergelen")
        fastq_minlen = self.option("fastq_minlen")
        fastq_maxee = self.option("fastq_maxee")
        fastq_truncqual = self.option("fastq_truncqual")
        fastq_maxns = self.option("fastq_maxns")
        fastq_ascii = self.option("fastq_ascii")
        threads = self.option("threads")
        #
        # Go
        # Opening the file to get all the lines stocked in the list csv_content
        fastq_and_fasta_list = [] # unique pairs of fastq files with the corresponding fasta file
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
                    fout_name = sample_info[7] + "_" + sample_info[4] + "_" +sample_info[6] + ".fasta"
                    PathManager.mkdir_p(fasta_dir)
                    fastainfo_line = line.strip() + '\t' + fout_name + '\n'
                    fastainfo_fout.write(fastainfo_line) # write fastainfo line
                    fasta_abspath = os.path.join(fasta_dir, fout_name)
                    if not (fastq_fw_abspath, fastq_rv_abspath, fasta_abspath) in fastq_and_fasta_list:
                        fastq_and_fasta_list.append((fastq_fw_abspath, fastq_rv_abspath, fasta_abspath))

        #########################################
        #
        # Loop of fastq pairs and run vsearch
        #
        #########################################
        for fastq_fw_abspath, fastq_rv_abspath, fasta_abspath in fastq_and_fasta_list:
            command = "vsearch -fastq_mergepairs {} --reverse {}".format(fastq_fw_abspath, fastq_rv_abspath)
            command += " --fastaout {} --fastq_minovlen {}".format(fasta_abspath, fastq_minovlen)
            command += " --fastq_maxmergelen {}".format(fastq_maxmergelen)
            command += " --fastq_minmergelen {}".format(fastq_minmergelen)
            command += " --fastq_minlen {}".format(fastq_minlen)
            command += " --fastq_maxee {}".format(fastq_maxee)
            command += " --fastq_truncqual {}".format(fastq_truncqual)
            command += " --fastq_maxns {}".format(fastq_maxns)
            command += " --fastq_ascii {}".format(fastq_ascii)
            command += " --threads {}".format(threads)
            Logger.instance().info(command)
            run_result = subprocess.run(command.split(), stdout=subprocess.PIPE)
            Logger.instance().info(run_result.stdout)
            try:
                assert run_result.returncode == 0 # vsearch exited correctly?
            except Exception:
                msg = "Vsearch error: The program has exited due to an error." \
                      "Verify the command and the input and output files."
                Logger.instance().error(run_result.stderr)
                Logger.instance().error(VTAMexception(msg))
                sys.exit(1)
