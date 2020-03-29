import multiprocessing
import os
import pandas
import pathlib
import sys
import yaml

from vtam.utils.Logger import Logger
from vtam.utils.SortReadsRunner import SortReadsRunner
from vtam.utils.VSearch import VSearch
from vtam.utils.VTAMexception import VTAMexception

from vtam.utils.constants import parameters_numerical

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


class CommandSortReads(object):
    """Class for the Merge command"""

    @classmethod
    def main(cls, fastainfo, fastadir, outdir, params=None, num_threads=multiprocessing.cpu_count()):

        min_id = parameters_numerical['min_id']
        minseqlength = parameters_numerical['minseqlength']
        overhang = parameters_numerical['overhang']

        if not (params is None):
            min_id = params['min_id']
            minseqlength = params['minseqlength']
            overhang = params['overhang']

        fastainfo_df = pandas.read_csv(fastainfo, sep='\t', header=0)

        pathlib.Path(outdir).mkdir(exist_ok=True)

        sorted_read_info_df = pandas.DataFrame()

        ################################################################################################################
        #
        # Loop over fasta files to sort reads per fasta
        #
        ################################################################################################################

        for fasta_filename in sorted(fastainfo_df.Fasta.unique().tolist()):

            Logger.instance().debug("Analysing FASTA file: {}".format(fasta_filename))

            fasta_path = os.path.join(fastadir, fasta_filename)
            fasta_info_df_i = fastainfo_df.loc[fastainfo_df.Fasta == fasta_filename]

            alignement_parameters = {'min_id': min_id, 'minseqlength': minseqlength, 'overhang': overhang}

            # Create SortReadsRunner
            sort_reads_runner = SortReadsRunner(fasta_information_df=fasta_info_df_i,
                                                fasta_path=fasta_path,
                                                alignement_parameters=alignement_parameters, outdir=outdir, num_threads=num_threads)

            sorted_read_info_df_i, sorted_read_dic = sort_reads_runner.run()
            sorted_read_info_df = sorted_read_info_df.append(sorted_read_info_df_i)

            ################################################################################################################
            #
            # Write to trimmed fasta
            #
            ################################################################################################################

            for sorted_read_filename in sorted_read_dic:

                sorted_read_list = sorted_read_dic[sorted_read_filename]
                sorted_read_path = os.path.join(outdir, sorted_read_filename)
                with open(sorted_read_path, "w") as fout:
                    fout.write("\n".join(sorted_read_list))

        fasta_trimmed_info_tsv = os.path.join(outdir, 'readinfo.tsv')
        sorted_read_info_df.to_csv(fasta_trimmed_info_tsv, sep="\t", header=True, index=False)
