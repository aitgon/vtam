import multiprocessing
import os
import pathlib
import sys
import yaml

from vtam.utils.VSearch import VSearch
from vtam.utils.VTAMexception import VTAMexception
from vtam.utils.Logger import Logger


class VSearchMergeRunner(object):

    parameters = {
    'fastq_ascii': 33,
    'fastq_maxee': 1,
    'fastq_maxmergelen': 500,
    'fastq_maxns': 0,
    'fastq_minlen': 50,
    'fastq_minmergelen': 100,
    'fastq_minovlen': 50,
    'fastq_truncqual': 10,
    "threads": int(multiprocessing.cpu_count()),
    }

    def __init__(self, fastq_fw_abspath, fastq_rv_abspath, fasta_abspath, params_yml=None, threads=None):

        self.fastq_fw_abspath = fastq_fw_abspath
        self.fastq_rv_abspath = fastq_rv_abspath
        self.fasta_abspath = fasta_abspath
        self.params_yml = params_yml
        self.threads = threads

    def load_parameters(self):

        self.parameters['fastq_mergepairs'] = self.fastq_fw_abspath
        self.parameters['reverse'] = self.fastq_rv_abspath
        self.parameters['fastaout'] = self.fasta_abspath
        if not (self.threads is None):
            self.parameters['threads'] = self.threads

        # Read parameters
        if not (self.params_yml is None):
            with open(self.params_yml, 'r') as fin:
                user_params_dic = yaml.load(fin, Loader=yaml.SafeLoader)
                for k_in in user_params_dic:
                    if k_in in VSearchMergeRunner.parameters.keys():
                        self.parameters[k_in] = user_params_dic[k_in]
                    else:
                        Logger.instance().error(
                            VTAMexception("One of the Merge parameters is not used by VSearch. "
                                          "Please set one or more of these parameters: {}".format(self.parameters.keys())))
                        sys.exit(1)

    def run(self):

        self.load_parameters()  # update vsearch parameters from user

        # Add double dash '--' to all parameters to create vsearch parameters
        vsearch_parameters = {}
        for par in self.parameters:
            vsearch_parameters['--{}'.format(par)] = self.parameters[par]
        vsearch_cluster = VSearch(parameters=vsearch_parameters)
        vsearch_cluster.run()


class CommandMerge(object):
    """Class for the Merge command"""

    @classmethod
    def main(cls, fastqinfo, fastqdir, fastainfo, fastadir, params=None, num_threads=multiprocessing.cpu_count()):
        #
        # Go
        # Opening the file to get all the lines stocked in the list csv_content
        fastq_and_fasta_list = [] # unique pairs of fastq files with the corresponding fasta_path file
        with open(fastqinfo, 'r') as csv_file:
            with open(fastainfo, 'w') as fastainfo_fout:
                next(csv_file) # skip header of fastqinfo
                fastainfo_header = "TagFwd	PrimerFwd	TagRev	PrimerRev	Marker	 " \
                                   "Biosample	Replicate	Run	FastqFwd	FastqRev	Fasta\n"
                fastainfo_fout.write(fastainfo_header) # write header of fastainfo
                for line in csv_file:
                    sample_info = line.strip().split("\t")
                    fastq_fw_abspath = os.path.join(fastqdir, sample_info[8])
                    fastq_rv_abspath = os.path.join(fastqdir, sample_info[9])
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
                    pathlib.Path(fastadir).mkdir(exist_ok=True)
                    fastainfo_line = line.strip() + '\t' + fasta_merged_basename + '\n'
                    fastainfo_fout.write(fastainfo_line) # write fastainfo line
                    fasta_abspath = os.path.join(fastadir, fasta_merged_basename)
                    if not (fastq_fw_abspath, fastq_rv_abspath, fasta_abspath) in fastq_and_fasta_list:
                        fastq_and_fasta_list.append((fastq_fw_abspath, fastq_rv_abspath, fasta_abspath))

        #########################################
        #
        # Loop of fastq pairs and run vsearch
        #
        #########################################
        for fastq_fw_abspath, fastq_rv_abspath, fasta_abspath in fastq_and_fasta_list:

            vsearch_merge_runner = VSearchMergeRunner(fastq_fw_abspath, fastq_rv_abspath, fasta_abspath,
                                                      params_yml=params, threads=num_threads)
            vsearch_merge_runner.load_parameters()
            vsearch_merge_runner.run()

