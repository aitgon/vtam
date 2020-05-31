import multiprocessing
import os
import pandas
import pathlib
import sys

from vtam.utils.ParamsFile import ParamsFile
from vtam.utils.SampleInformationFile import SampleInformationFile
from vtam.utils.VSearch import VSearch
from vtam.utils.VTAMexception import VTAMexception
from vtam.utils.Logger import Logger


class CommandMerge(object):
    """Class for the Merge command"""

    @staticmethod
    def main(fastqinfo, fastqdir, fastainfo, fastadir, params=None, num_threads=multiprocessing.cpu_count()):

        ############################################################################################
        #
        # params.yml parameters
        #
        ############################################################################################

        params_dic = ParamsFile(params).get_params_dic()

        ############################################################################################
        #
        # Read fastq information into df
        #
        ############################################################################################

        fastqinfo_df = SampleInformationFile(fastqinfo).read_tsv_into_df()

        pathlib.Path(
            os.path.dirname(fastainfo)).mkdir(
            parents=True,
            exist_ok=True)
        pathlib.Path(fastadir).mkdir(parents=True, exist_ok=True)

        fastainfo_df = pandas.DataFrame()

        ############################################################################################
        #
        # Loop over fastq pairs to merge
        #
        ############################################################################################

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

            ########################################################################################
            #
            # Run vsearch merge
            #
            ########################################################################################

            vsearch_args_dic = {}

            vsearch_args_dic['fastq_ascii'] = params_dic['fastq_ascii']
            vsearch_args_dic['fastq_maxee'] = params_dic['fastq_maxee']
            vsearch_args_dic['fastq_maxmergelen'] = params_dic['fastq_maxmergelen']
            vsearch_args_dic['fastq_maxns'] = params_dic['fastq_maxns']
            vsearch_args_dic['fastq_minlen'] = params_dic['fastq_minlen']
            vsearch_args_dic['fastq_minmergelen'] = params_dic['fastq_minmergelen']
            vsearch_args_dic['fastq_minovlen'] = params_dic['fastq_minovlen']
            vsearch_args_dic['fastq_truncqual'] = params_dic['fastq_truncqual']

            vsearch_args_dic['fastq_mergepairs'] = fastq_fw_abspath
            vsearch_args_dic['reverse'] = fastq_rv_abspath
            vsearch_args_dic['fastaout'] = out_fasta_path
            vsearch_args_dic['threads'] = num_threads

            vsearch_cluster = VSearch(parameters=vsearch_args_dic)
            vsearch_cluster.run()

            fastq_info_df_i = fastq_info_df_i[['run', 'marker', 'biosample', 'replicate', 'tagfwd',
                                               'primerfwd', 'tagrev', 'primerrev']]
            fastq_info_df_i['mergedfasta'] = fasta_merged_basename
            fastainfo_df = pandas.concat(
                [fastainfo_df, fastq_info_df_i], axis=0)

        fastainfo_df.to_csv(fastainfo, sep="\t", header=True, index=False)
