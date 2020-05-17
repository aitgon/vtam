import multiprocessing
import os
import pandas
import pathlib
import sys

from vtam.utils.SampleInformationFile import SampleInformationFile
from vtam.utils.VSearch import VSearch
from vtam.utils.VTAMexception import VTAMexception
from vtam.utils.Logger import Logger
from vtam.utils import constants


class CommandMerge(object):
    """Class for the Merge command"""

    @staticmethod
    def main(
            fastqinfo,
            fastqdir,
            fastainfo,
            fastadir,
            params=None,
            num_threads=multiprocessing.cpu_count()):

        #######################################################################
        #
        # Parameters
        #
        #######################################################################

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
                    merge_vsearch_parameters[parameter_i] = params[parameter_i]

        #######################################################################
        #
        # Read fastq information into df
        #
        #######################################################################

        fastqinfo_df = SampleInformationFile(fastqinfo).read_tsv_into_df()

        pathlib.Path(
            os.path.dirname(fastainfo)).mkdir(
            parents=True,
            exist_ok=True)
        pathlib.Path(fastadir).mkdir(parents=True, exist_ok=True)

        fastainfo_df = pandas.DataFrame()

        #######################################################################
        #
        # Loop over fastq pairs to merge
        #
        #######################################################################

        for fastqfwd, fastqrev in fastqinfo_df[[
                'fastqfwd', 'fastqrev']].drop_duplicates().values:

            # fastq_info_series = fastqinfo_df.loc[]
            # fastq_info_df_i = fastq_info_series.to_frame().T
            fastq_info_df_i = fastqinfo_df.loc[(fastqinfo_df.fastqfwd == fastqfwd) & (
                fastqinfo_df.fastqrev == fastqrev)]

            fastq_fw_abspath = os.path.join(fastqdir, fastqfwd)
            fastq_rv_abspath = os.path.join(fastqdir, fastqrev)

            Logger.instance().debug(
                "Analysing FASTQ files: {} and ".format(
                    fastqfwd, fastqrev))

            try:
                pathlib.Path(fastq_fw_abspath).resolve(strict=True)
            except FileNotFoundError:
                Logger.instance().error(
                    VTAMexception(
                        "VTAMexception: This FASTQ file was not found: {}.".format(fastq_fw_abspath)))
                sys.exit(1)
            try:
                pathlib.Path(fastq_rv_abspath).resolve(strict=True)
            except FileNotFoundError:
                Logger.instance().error(
                    VTAMexception(
                        "VTAMexception: This FASTQ file was not found: {}.".format(fastq_rv_abspath)))
                sys.exit(1)

            # fasta_merged_basename = '.'.join(os.tsv_path.basename(fastq_fw_abspath).split('.')[0:-1]) + '_merged.fasta'
            fasta_merged_basename = os.path.basename(
                fastq_fw_abspath).replace('.fastq', '.fasta')
            out_fasta_path = os.path.join(fastadir, fasta_merged_basename)

            ###################################################################
            #
            # Run vsearch merge
            #
            ###################################################################

            # vsearch_parameters = {}
            # for par in merge_vsearch_parameters:
            #     vsearch_parameters['--{}'.format(par)] = self.parameters[par]
            merge_vsearch_parameters['fastq_mergepairs'] = fastq_fw_abspath
            merge_vsearch_parameters['reverse'] = fastq_rv_abspath
            merge_vsearch_parameters['fastaout'] = out_fasta_path
            merge_vsearch_parameters['threads'] = num_threads
            vsearch_cluster = VSearch(parameters=merge_vsearch_parameters)
            vsearch_cluster.run()

            fastq_info_df_i = fastq_info_df_i[['run',
                                               'marker',
                                               'biosample',
                                               'replicate',
                                               'tagfwd',
                                               'primerfwd',
                                               'tagrev',
                                               'primerrev']]
            fastq_info_df_i['mergedfasta'] = fasta_merged_basename
            fastainfo_df = pandas.concat(
                [fastainfo_df, fastq_info_df_i], axis=0)

        fastainfo_df.to_csv(fastainfo, sep="\t", header=True, index=False)
