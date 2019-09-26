#!/usr/bin/env python

import argparse
import os

from wopmetabarcoding.utils.PathFinder import PathFinder

class ArgParser():

    @staticmethod
    def get_arg_parser_base():
        """These are base VTAM arguments"""
        parser = argparse.ArgumentParser(
                    description='vtam',
                    usage='''vtam <command> [<args>]

        The available vtam commands are:
           merge     Merge FASTQ files
           otu       Sort reads, creates variants, filters variants and creates OTU table
           optimize  Shows positive, negative and other variants frequencies to help select parameters 
        ''')
        ###################
        #
        # required
        #
        ###################
        parser.add_argument('command', help='Subcommand to run: merge or otu or optimize',
                            choices=['merge', 'otu', 'optimize'])
        ###################
        #
        # optional
        #
        ###################
        parser.add_argument('--dry-run', '-n', dest='dryrun', action='store_true',
                            help="Only display what would have been done.", required=False)
        parser.add_argument('-F', '--forceall', dest='forceall', action='store_true', help="Force argument of WopMars",
                            required=False)
        parser.add_argument('--log', nargs=1, dest='log',
                            help="Write log to file.", required=False,
                            type=os.path.abspath)
        parser.add_argument('--params', nargs=1, help="YML file with parameter values", required=False)
        parser.add_argument('-t', '--targetrule', nargs=1, help="Execute the workflow to the given target "
                                                                "RULE: SampleInformation, ...", required=False)
        parser.add_argument('-v', dest='verbose', action='count',
                            help="Set verbosity level, eg. None (Error level) -v (Info level) or -vv (Debug level)",
                            required=False)
        return parser


    @staticmethod
    def get_arg_parser_merge():
        """These are 'vtam merge' arguments"""
        # Base parser
        parser = ArgParser.get_arg_parser_base()
        ###################
        #
        # required
        #
        ###################
        parser.add_argument('--db', dest='db', nargs=1, help="SQLITE file with DB", required=True,
                            type=os.path.abspath)
        parser.add_argument('--fastqinfo', dest='fastqinfo', nargs=1, help="TSV file with FASTQ sample information",
                            required=True, type=lambda x:
                            PathFinder.check_file_exists_and_is_nonempty(x, "Verify the '--fastqinfo' argument"))
        parser.add_argument('--fastainfo', dest='fastainfo', nargs=1, help="TSV file with FASTA sample information",
                            required=True, type=os.path.abspath)
        parser.add_argument('--fastqdir', dest='fastqdir', nargs=1, help="Directory with FASTQ files", required=True,
                            type=lambda x:
                            PathFinder.check_dir_exists_and_is_nonempty(x, "Verify the '--fastqdir' argument"))
        parser.add_argument('--fastadir', dest='fastadir', nargs=1, help="Directory with FASTA files", required=True,
                            type=os.path.abspath)
        return parser


    @staticmethod
    def get_arg_parser_otu():
        """These are 'vtam merge' arguments"""
        # Base parser
        parser = ArgParser.get_arg_parser_base()
        ###################
        #
        # required
        #
        ###################
        parser.add_argument('--db', dest='db', nargs=1, help="SQLITE file with DB", required=True,
                            type=os.path.abspath)
        parser.add_argument('--fastainfo', dest='fastainfo', nargs=1, help="TSV file with FASTA sample information",
                            type=lambda x:
                            PathFinder.check_file_exists_and_is_nonempty(x, "Verify the '--fastainfo' argument")
                            , required=True)
        parser.add_argument('--fastadir', dest='fastadir', nargs=1, help="Directory with FASTA files",
                            type=lambda x:
                            PathFinder.check_dir_exists_and_is_nonempty(x, "Verify the '--fastadir' argument"),
                            required=True)
        parser.add_argument('--outdir', nargs=1, help="Directory for output", required=True)
        # parser.add_argument('--filter_lfn_variant', nargs=1,
        #                     help="Boolean 0|1 to filter_lfn_variant (1) or filter_lfn_variant_replicate (0)",
        #                     required=True)
        parser.add_argument('--filter_lfn_variant', default=False, action='store_true',
            help="If present, VTAM runs filter_lfn_variant. Otherwise, VTAM runs filter_lfn_variant_replicate.")
        parser.add_argument('--threshold_specific', nargs=1, help="Variant or variant-replicate specific threshold", required=False)
        return parser

