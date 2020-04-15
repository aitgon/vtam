import multiprocessing
import os
import pandas
import pathlib
import sys
import yaml

from vtam.utils.PathManager import PathManager
from vtam.utils.VSearch import VSearch
from vtam.utils.VTAMexception import VTAMexception
from vtam.utils.Logger import Logger
from vtam.utils import constants


# class VSearchMergeRunner(object):
#
#     # parameters = {
#     # 'fastq_ascii': 33,
#     # 'fastq_maxee': 1,
#     # 'fastq_maxmergelen': 500,
#     # 'fastq_maxns': 0,
#     # 'fastq_minlen': 50,
#     # 'fastq_minmergelen': 100,
#     # 'fastq_minovlen': 50,
#     # 'fastq_truncqual': 10,
#     # "threads": int(multiprocessing.cpu_count()),
#     # }
#
#     def __init__(self, fastq_fw_abspath, fastq_rv_abspath, fasta_abspath, params_yml=None, threads=None):
#
#         self.fastq_fw_abspath = fastq_fw_abspath
#         self.fastq_rv_abspath = fastq_rv_abspath
#         self.fasta_abspath = fasta_abspath
#         self.params_yml = params_yml
#         self.threads = threads
#
#     def load_parameters(self):
#
#         self.parameters['fastq_mergepairs'] = self.fastq_fw_abspath
#         self.parameters['reverse'] = self.fastq_rv_abspath
#         self.parameters['fastaout'] = self.fasta_abspath
#         if not (self.threads is None):
#             self.parameters['threads'] = self.threads
#
#         # Read parameters
#         if not (self.params_yml is None):
#             with open(self.params_yml, 'r') as fin:
#                 user_params_dic = yaml.load(fin, Loader=yaml.SafeLoader)
#                 for k_in in user_params_dic:
#                     if k_in in VSearchMergeRunner.parameters.keys():
#                         self.parameters[k_in] = user_params_dic[k_in]
#                     else:
#                         Logger.instance().error(
#                             VTAMexception("One of the Merge parameters is not used by VSearch. "
#                                           "Please set one or more of these parameters: {}".format(self.parameters.keys())))
#                         sys.exit(1)
#
#     def run(self):
#
#         self.load_parameters()  # update vsearch parameters from user
#
#         # Add double dash '--' to all parameters to create vsearch parameters
#         vsearch_parameters = {}
#         for par in self.parameters:
#             vsearch_parameters['--{}'.format(par)] = self.parameters[par]
#         vsearch_cluster = VSearch(parameters=vsearch_parameters)
#         vsearch_cluster.run()


class CommandMerge(object):
    """Class for the Merge command"""

    @classmethod
    def main(cls, fastqinfo, fastqdir, fastainfo, fastadir, params=None, num_threads=multiprocessing.cpu_count()):

        ################################################################################################################
        #
        # Parameters
        #
        ################################################################################################################

        merge_vsearch_parameters = {}

        params_default = constants.get_dic_params_default()
        merge_vsearch_parameters['fastq_ascii'] = params_default['fastq_ascii']
        merge_vsearch_parameters['fastq_maxee'] = params_default['fastq_maxee']
        merge_vsearch_parameters['fastq_maxmergelen'] = params_default['fastq_maxmergelen']
        merge_vsearch_parameters['fastq_maxns'] = params_default['fastq_maxns']
        merge_vsearch_parameters['fastq_minlen'] = params_default['fastq_minlen']
        merge_vsearch_parameters['fastq_minmergelen'] = params_default['fastq_minmergelen']
        merge_vsearch_parameters['fastq_minovlen'] = params_default['fastq_minovlen']
        merge_vsearch_parameters['fastq_truncqual'] = params_default['fastq_truncqual']

        if not (params is None):
            for parameter_i in merge_vsearch_parameters:
                if parameter_i in params:
                    merge_vsearch_parameters['vsearch_fastq_ascii'] = params[parameter_i]

        ################################################################################################################
        #
        # Open fastq information
        #
        ################################################################################################################

        fastqinfo_df = pandas.read_csv(fastqinfo, sep='\t', header=0)
        fastqinfo_df.columns = fastqinfo_df.columns.str.lower()

        pathlib.Path(fastadir).mkdir(parents=True, exist_ok=True)

        # tempdir = PathManager.instance().get_tempdir()
        fastainfo_df = pandas.DataFrame()

        ################################################################################################################
        #
        # Loop over fastq pairs to merge
        #
        ################################################################################################################

        for fastqfwd, fastqrev in fastqinfo_df[['fastqfwd', 'fastqrev']].drop_duplicates().values:

            # fastq_info_series = fastqinfo_df.loc[]
            # fastq_info_df_i = fastq_info_series.to_frame().T
            fastq_info_df_i = fastqinfo_df.loc[(fastqinfo_df.fastqfwd == fastqfwd) & (fastqinfo_df.fastqrev == fastqrev)]

            fastq_fw_abspath = os.path.join(fastqdir, fastqfwd)
            fastq_rv_abspath = os.path.join(fastqdir, fastqrev)

            Logger.instance().debug("Analysing FASTQ files: {} and ".format(fastqfwd, fastqrev))

            try:
                pathlib.Path(fastq_fw_abspath).resolve(strict=True)
            except FileNotFoundError:
                Logger.instance().error(
                    VTAMexception("VTAMexception: This FASTQ file was not found: {}.".format(fastq_fw_abspath)))
                sys.exit(1)
            try:
                pathlib.Path(fastq_rv_abspath).resolve(strict=True)
            except FileNotFoundError:
                Logger.instance().error(
                    VTAMexception("VTAMexception: This FASTQ file was not found: {}.".format(fastq_rv_abspath)))
                sys.exit(1)

            # fasta_merged_basename = '.'.join(os.path.basename(fastq_fw_abspath).split('.')[0:-1]) + '_merged.fasta'
            fasta_merged_basename = os.path.basename(fastq_fw_abspath).replace('.fastq', '.fasta')
            out_fasta_path = os.path.join(fastadir, fasta_merged_basename)

            ################################################################################################################
            #
            # Run vsearch merge
            #
            ################################################################################################################

            # vsearch_parameters = {}
            # for par in merge_vsearch_parameters:
            #     vsearch_parameters['--{}'.format(par)] = self.parameters[par]
            merge_vsearch_parameters['fastq_mergepairs'] = fastq_fw_abspath
            merge_vsearch_parameters['reverse'] = fastq_rv_abspath
            merge_vsearch_parameters['fastaout'] = out_fasta_path
            merge_vsearch_parameters['threads'] = num_threads
            vsearch_cluster = VSearch(parameters=merge_vsearch_parameters)
            vsearch_cluster.run()

            fastq_info_df_i = fastq_info_df_i[['run', 'marker', 'biosample', 'replicate', 'tagfwd', 'primerfwd', 'tagrev', 'primerrev']]
            fastq_info_df_i['mergedfasta'] = fasta_merged_basename
            fastainfo_df = pandas.concat([fastainfo_df, fastq_info_df_i], axis=0)

        fastainfo_df.to_csv(fastainfo, sep="\t", header=True, index=False)

