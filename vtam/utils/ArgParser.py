#!/usr/bin/env python

import argparse
import os

from vtam.utils.Logger import Logger
from vtam.utils.OptionManager import OptionManager
from vtam.utils.VTAMexception import VTAMexception
from vtam.utils.PathManager import PathManager


class ArgParser():

    @staticmethod
    def verify_an_store_blast_db_argument(blast_db_arg, is_abspath=False):
        """Verifies --blast_db argument. Must be exactly two arguments.

        - First argument: REQUIRED. Blast DB directory
        - Second argument: REQUIRED. Basename of indexed Blast files that is 'nt' for full Blast DB. If basename not 'nt',
        - Third argument: Optional. If custom Blast DB, the this argument must give the path to a TSV file with the
        mapping from sequence IDs to tax IDs.

        :param error_message: Optional message to help debug the problem
        :param abspath: If True, returns abspath
        :return: void
        """
        if not 'blast_db_dir' in OptionManager.instance(): # First argument of blast_db, which is the blast db dir
            if is_abspath:
                blast_db_arg = os.path.abspath(blast_db_arg)
            PathManager.check_dir_exists_and_is_nonempty(blast_db_arg,
                                                         error_message='Verify the existance and non-emptyness of the directory in '
                                                                       'the first argument to the --blast_db argument',
                                                         is_abspath=is_abspath)
            OptionManager.instance()['blast_db_dir'] = blast_db_arg
        else:
            OptionManager.instance()['blast_db_basename'] = blast_db_arg
            ###############################
            #
            # blast_nhr
            #
            ###############################
            # blast_nhr = os.path.join(OptionManager.instance()['blast_db_dir'], blast_db_arg + '.nhr')
            # Verify presence of at least one copy of these files: nhr, nin, nog, nsd, nsi, nsq
            one_file_exists = {'nhr': 0, 'nin': 0, 'nog': 0, 'nsd': 0, 'nsi': 0, 'nsq': 0}
            for fname in os.listdir(OptionManager.instance()['blast_db_dir']):
                if fname.endswith('.nhr'):
                    one_file_exists['nhr'] = 1
                elif fname.endswith('.nin'):
                    one_file_exists['nin'] = 1
                elif fname.endswith('.nog'):
                    one_file_exists['nog'] = 1
                elif fname.endswith('.nsd'):
                    one_file_exists['nsd'] = 1
                elif fname.endswith('.nsi'):
                    one_file_exists['nsi'] = 1
                elif fname.endswith('.nsq'):
                    one_file_exists['nsq'] = 1
            if not sum(one_file_exists.values()) == 6:
                raise Logger.instance().error(VTAMexception("Verify if there are NHR, NIN, NOG, NSD, NSI and NSQ files "
                                                            "in the blast directory"))
            # # do stuff if a file .true doesn't exist.
            # PathManager.check_file_exists_and_is_nonempty(blast_nhr,
            #                                              error_message='Verify the existance of Blast files defined by the basename '
            #                                                            'in the second argument to the --blast_db argument',
            #                                              is_abspath=is_abspath)
            # OptionManager.instance()['blast_nhr'] = blast_nhr
            # ###############################
            # #
            # # blast_nin
            # #
            # ###############################
            # # blast_nin = os.path.join(OptionManager.instance()['blast_db_dir'], blast_db_arg + '.nin')
            # PathManager.check_file_exists_and_is_nonempty(blast_nin,
            #                                              error_message='Verify the existance of Blast files defined by the basename '
            #                                                            'in the second argument to the --blast_db argument',
            #                                              is_abspath=is_abspath)
            # # OptionManager.instance()['blast_nin'] = blast_nin
            # ###############################
            # #
            # # blast_nog

            # #
            # ###############################
            # # blast_nog = os.path.join(OptionManager.instance()['blast_db_dir'], blast_db_arg + '.nog')
            # PathManager.check_file_exists_and_is_nonempty(blast_nog,
            #                                              error_message='Verify the existance of Blast files defined by the basename '
            #                                                            'in the second argument to the --blast_db argument',
            #                                              is_abspath=is_abspath)
            # # OptionManager.instance()['blast_nog'] = blast_nog
            # ###############################
            # #
            # # blast_nsd
            # #
            # ###############################
            # # blast_nsd = os.path.join(OptionManager.instance()['blast_db_dir'], blast_db_arg + '.nsd')
            # PathManager.check_file_exists_and_is_nonempty(blast_nsd,
            #                                              error_message='Verify the existance of Blast files defined by the basename '
            #                                                            'in the second argument to the --blast_db argument',
            #                                              is_abspath=is_abspath)
            # # OptionManager.instance()['blast_nsd'] = blast_nsd
            # ###############################
            # #
            # # blast_nsi
            # #
            # ###############################
            # # blast_nsi = os.path.join(OptionManager.instance()['blast_db_dir'], blast_db_arg + '.nsi')
            # blast_nsi = PathManager.check_file_exists_and_is_nonempty(blast_nsi,
            #                                              error_message='Verify the existance of Blast files defined by the basename '
            #                                                            'in the second argument to the --blast_db argument',
            #                                              is_abspath=is_abspath)
            # # OptionManager.instance()['blast_nsi'] = blast_nsi
            # ###############################
            # #
            # # blast_nsq
            # #
            # ###############################
            # # blast_nsq = os.path.join(OptionManager.instance()['blast_db_dir'], blast_db_arg + '.nsq')
            # PathManager.check_file_exists_and_is_nonempty(blast_nsq,
            #                                              error_message='Verify the existance of Blast files defined by the basename '
            #                                                            'in the second argument to the --blast_db argument',
            #                                              is_abspath=is_abspath)
            # # OptionManager.instance()['blast_nsq'] = blast_nsq
        return blast_db_arg

    @staticmethod
    def get_arg_parser(is_abspath=False):
        """

        :param is_abspath: If True, return absolute paths
        :return:
        """
        # create the top-level parser
        parser_vtam = argparse.ArgumentParser(add_help=False)
        parser_vtam.add_argument('--db', action='store', default='db.sqlite', required=False,
                                 type=lambda x: os.path.abspath(x) if is_abspath else x, help="SQLITE file with DB")
        parser_vtam.add_argument('--dry-run', '-n', dest='dryrun', action='store_true', required=False,
                                 help="Only display what would have been done.")
        parser_vtam.add_argument('-F', '--forceall', dest='forceall', action='store_true',
                                 help="Force argument of WopMars", required=False)
        parser_vtam.add_argument('--log', dest='log_file', action='store', help="Write log to file.", required=False,
                                 type=lambda x: os.path.abspath(x) if is_abspath else x)
        parser_vtam.add_argument('--params', action='store', default=None, help="YML file with parameter values",
                                 required=False)
        parser_vtam.add_argument('-t', '--targetrule', dest='targetrule', action='store', default=None,
                                 help="Execute the workflow to the given target RULE: SampleInformation, ...",
                                 required=False)
        parser_vtam.add_argument('-f', '--sourcerule', dest='sourcerule', action='store', default=None,
                                 help="Execute the workflow from the given RULE.",
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
                                                                                                    error_message="Verify the '--fastqinfo' argument",
                                                                                                    is_abspath=is_abspath))
        parser_vtam_merge.add_argument('--fastainfo', action='store', help="TSV file with FASTA sample information",
                                       required=True, type=lambda x: os.path.abspath(x) if is_abspath else x)
        parser_vtam_merge.add_argument('--fastqdir', action='store', help="Directory with FASTQ files", required=True,
                                       type=lambda x:
                                       PathManager.check_dir_exists_and_is_nonempty(x,
                                                                                    error_message="Verify the '--fastqdir' argument",
                                                                                    is_abspath=is_abspath))
        parser_vtam_merge.add_argument('--fastadir', action='store', help="Directory with FASTA files", required=True,
                                       type=lambda x: os.path.abspath(x) if is_abspath else x)
        parser_vtam_merge.set_defaults(command='merge')  # This attribute will trigget the good command

        #############################################
        #
        # create the parser for the "otu" command
        #
        #############################################
        parser_vtam_otu = subparsers.add_parser('otu', add_help=True, parents=[parser_vtam])
        parser_vtam_otu.add_argument('--fastainfo', action='store',
                                     help="REQUIRED: TSV file with FASTA sample information",
                                     required=True, type=lambda x:
            PathManager.check_file_exists_and_is_nonempty(x,
                                                          error_message="Verify the '--fastainfo' argument",
                                                          is_abspath=is_abspath))
        parser_vtam_otu.add_argument('--fastadir', action='store', help="REQUIRED: Directory with FASTA files",
                                     required=True,
                                     type=lambda x:
                                     PathManager.check_file_exists_and_is_nonempty(x,
                                                                                   error_message="Verify the '--fastadir' argument",
                                                                                   is_abspath=is_abspath))
        parser_vtam_otu.add_argument('--outdir', action='store', help="REQUIRED: Directory for output", default="out",
                                     required=True)
        parser_vtam_otu.add_argument('--blast_db', action='store', nargs=2,
                                     help="REQUIRED: One argument, optional the second argument"
                                          "-First argument: Blast DB directory (Full or custom one)"
                                          "-Second argument: Blast DB file basename", default="blast_dir",
                                     required=True,
                                     type=lambda x:
                                     ArgParser.verify_an_store_blast_db_argument(x, is_abspath=is_abspath))
        parser_vtam_otu.add_argument('--map_taxids', action='store',
                                     help="TSV file with mapping from NCBI sequence IDs to tax IDs."
                                          "Required if working with custome DB.",
                                     required=False, type=lambda x:
            PathManager.check_file_exists_and_is_nonempty(x,
                                                          error_message="Verify the file in the '--map_taxids' argument",
                                                          is_abspath=is_abspath))
        parser_vtam_otu.add_argument('--taxonomy', dest='taxonomy', action='store',
                                     help="""REQUIRED: SQLITE DB with taxonomy information.

        This database is create with the command: create_db_taxonomy. For instance

        create_db_taxonomy -o taxonomy.sqlite to create a database in the current directory.""",
                                     required=True,
                                     type=lambda x: PathManager.check_file_exists_and_is_nonempty(x,
                                                                  error_message="Verify the '--taxonomy' argument",
                                                                  is_abspath=is_abspath))
        parser_vtam_otu.add_argument('--threshold_specific', default=None, action='store', required=False,
                                     help="TSV file with variant (col1: variant; col2: threshold) or variant-replicate "
                                          "(col1: variant; col2: replicate; col3: threshold)specific thresholds. Header expected.",
                                     type=lambda x: PathManager.check_file_exists_and_is_nonempty(x,
                                                                  error_message="Verify the '--threshold_specific' argument",
                                                                  is_abspath=is_abspath))
        parser_vtam_otu.set_defaults(command='otu')  # This attribute will trigget the good command
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
        parser_vtam_optimize.add_argument('--fastadir', action='store', help="REQUIRED:Directory with FASTA files",
                                          required=True,
                                          type=lambda x:
                                          PathManager.check_file_exists_and_is_nonempty(x,
                                                                                        error_message="Verify the '--fastadir' argument",
                                                                                        is_abspath=is_abspath))
        parser_vtam_optimize.add_argument('--outdir', action='store', help="Directory for output", default="out",
                                          required=True)
        parser_vtam_optimize.add_argument('--variant_known', action='store', help="TSV file with known variants",
                                          required=True)
        parser_vtam_optimize.set_defaults(command='optimize')  # This attribute will trigget the good command

        return parser_vtam
