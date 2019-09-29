#!/usr/bin/env python

import argparse
import os

from wopmetabarcoding.utils.PathManager import PathManager

class ArgParser():

    @staticmethod
    def get_arg_parser(abspath=False):
        """

        :param abspath: If True, return absolute paths
        :return:
        """
        # create the top-level parser
        parser_vtam = argparse.ArgumentParser(add_help=False)
        parser_vtam.add_argument('--db', action='store', default='db.sqlite', required=False,
                                 type=lambda x: os.path.abspath(x) if abspath else x, help="SQLITE file with DB")
        parser_vtam.add_argument('--dry-run', '-n', dest='dryrun', action='store_true', required=False,
                                 help="Only display what would have been done.")
        parser_vtam.add_argument('-F', '--forceall', dest='forceall', action='store_true',
                                 help="Force argument of WopMars", required=False)
        parser_vtam.add_argument('--log', dest='log_file', action='store', help="Write log to file.", required=False,
                                 type=lambda x: os.path.abspath(x) if abspath else x)
        parser_vtam.add_argument('--params', action='store', default=None, help="YML file with parameter values",
                                 required=False)
        parser_vtam.add_argument('-t', '--targetrule', dest='targetrule', action='store', default=None,
                                 help="Execute the workflow to the given target RULE: SampleInformation, ...",
                                 required=False)
        parser_vtam.add_argument('-v', dest='log_verbosity', action='count', default=0, required=False,
                                 help="Set verbosity level, eg. None (Error level) -v (Info level) or -vv (Debug level)"
                            )
        subparsers = parser_vtam.add_subparsers()

        #############################################
        #
        # create the parser for the "merge" command
        #
        #############################################
        parser_vtam_merge = subparsers.add_parser('merge', add_help=True, parents=[parser_vtam])
        parser_vtam_merge.add_argument('--fastqinfo', action='store', help="TSV file with FASTQ sample information",
                                       required=True,
                                       type=lambda x: PathManager.check_file_exists_and_is_nonempty(x,
                                                 error_message="Verify the '--fastqinfo' argument", abspath=abspath))
        parser_vtam_merge.add_argument('--fastainfo', action='store', help="TSV file with FASTA sample information",
                            required=True, type = lambda x: os.path.abspath(x) if abspath else x)
        parser_vtam_merge.add_argument('--fastqdir', action='store', help="Directory with FASTQ files", required=True,
                                       type=lambda x:
                                            PathManager.check_dir_exists_and_is_nonempty(x,
                                            error_message="Verify the '--fastqdir' argument", abspath=abspath))
        parser_vtam_merge.add_argument('--fastadir', action='store', help="Directory with FASTA files", required=True,
                                       type=lambda x: os.path.abspath(x) if abspath else x)
        parser_vtam_merge.set_defaults(command='merge') # This attribute will trigget the good command

        #############################################
        #
        # create the parser for the "otu" command
        #
        #############################################
        parser_vtam_otu = subparsers.add_parser('otu', add_help=True, parents=[parser_vtam])
        parser_vtam_otu.add_argument('--fastainfo', action='store', help="TSV file with FASTA sample information",
                                       required=True, type=lambda x:
                                            PathManager.check_file_exists_and_is_nonempty(x,
                                             error_message="Verify the '--fastainfo' argument", abspath=abspath))
        parser_vtam_otu.add_argument('--fastadir', action='store', help="Directory with FASTA files", required=True,
                                       type=lambda x:
                                            PathManager.check_file_exists_and_is_nonempty(x,
                                            error_message="Verify the '--fastadir' argument", abspath=abspath))
        parser_vtam_otu.add_argument('--outdir', action='store', help="Directory for output", default="out",
                                     required=False)
        parser_vtam_otu.add_argument('--filter_lfn_variant', default=False, action='store_true', required=False,
                    help="Boolean 0|1 to filter_lfn_variant (1) or filter_lfn_variant_replicate (0)")
        parser_vtam_otu.add_argument('--threshold_specific', default=False, action='store_true', required=False,
                                     help="Variant or variant-replicate specific threshold")
        parser_vtam_otu.set_defaults(command='otu') # This attribute will trigget the good command
        #############################################
        #
        # create the parser for the "optimize" command
        #
        #############################################
        parser_vtam_optimize = subparsers.add_parser('optimize', add_help=True, parents=[parser_vtam])
        parser_vtam_optimize.add_argument('--fastainfo', action='store', help="TSV file with FASTA sample information",
                                       required=True, type=lambda x:
                            PathManager.check_file_exists_and_is_nonempty(x,
                                                                 error_message="Verify the '--fastainfo' argument"))
        parser_vtam_optimize.add_argument('--fastadir', action='store', help="Directory with FASTA files", required=True,
                                       type=lambda x:
                                       PathManager.check_file_exists_and_is_nonempty(x,
                                                                    error_message="Verify the '--fastadir' argument"))
        parser_vtam_optimize.add_argument('--outdir', action='store', help="Directory for output", default="out",
                                     required=False)
        parser_vtam_optimize.add_argument('--variant_known', action='store', help="TSV file with known variants",
                                          required=True)
        parser_vtam_optimize.set_defaults(command='optimize') # This attribute will trigget the good command


        return parser_vtam

