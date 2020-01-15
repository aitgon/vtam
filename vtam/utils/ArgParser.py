import argparse
import multiprocessing
import os

from vtam.utils.Logger import Logger
from vtam.utils.OptionManager import OptionManager
from vtam.utils.VTAMexception import VTAMexception


class ArgParserChecker():
    """Methods to check arguments"""

    @staticmethod
    def check_parser_vtam_pool_markers_arg_runmarker(runmarker_tsv):
        # Check header and separator
        with open(runmarker_tsv, 'r') as fin:
            header_fields = fin.readline().strip().split('\t')
        if not set(header_fields) == {'run_name', 'marker_name'}:
            raise Logger.instance().error(VTAMexception("Verify the '--runmarker' argument"))
        return runmarker_tsv

    @staticmethod
    def check_file_exists_and_is_nonempty(path, error_message=None):
        """Checks if file exists and is not empty

        :param error_message: Optional message to help debug the problem
        :return: void
        """
        try:
            assert os.stat(path).st_size > 0
        except AssertionError as err:
            raise Logger.instance().error(VTAMexception("{}: {}".format(err, error_message)))
        except FileNotFoundError as err:
            raise Logger.instance().error(VTAMexception("{}: {}".format(err, error_message)))
        return path

    @staticmethod
    def check_dir_exists_and_is_nonempty(path, error_message=None):
        """Checks if directory exists and is not empty

        :param error_message: Optional message to help debug the problem
        :return: void
        """
        try:
            assert len(os.listdir(path)) > 0
            # assert True
        except AssertionError as err:
            raise Logger.instance().error(VTAMexception("{}: {}".format(err, error_message)))
        return path

    @staticmethod
    def check_blast_db_argument(blast_db):
        """Verifies --blast_db argument. Must be exactly two arguments.

        - First argument: REQUIRED. Blast DB directory

        :param error_message: Optional message to help debug the problem
        :param abspath: If True, returns abspath
        :return: void
        """

        ArgParserChecker.check_dir_exists_and_is_nonempty(blast_db,
                                             error_message='Verify the existance and non-emptyness of the directory in '
                                                           'the first argument to the --blast_db argument')
        OptionManager.instance()['blast_db'] = blast_db

        one_file_exists = {'nhr': 0, 'nin': 0, 'nog': 0, 'nsd': 0, 'nsi': 0, 'nsq': 0}
        for fname in os.listdir(OptionManager.instance()['blast_db']):
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
            raise Logger.instance().error(VTAMexception("Verify the Blast DB directory for files with 'nt' file name"
                                                        " and 'nhr', 'nin', 'nog', 'nsd', 'nsi' and 'nsq' file types."))
        return blast_db


class ArgParser:

    args_db = {'dest': 'db', 'action': 'store', 'default': 'db.sqlite', 'required': False,
               'help': "Database in SQLITE format"}
    args_log_file = {'dest': 'log_file', 'action': 'store', 'help': "Write log to file.", 'required': False}
    args_log_verbosity = {'dest': 'log_verbosity', 'action': 'count', 'default': 0, 'required': False,
                          'help': "Set verbosity level, eg. None (Error level) -v (Info level)."}

    @classmethod
    def get_arg_parser(cls):
        """

        :return:
        """
        # create the top-level parser
        parser_vtam = argparse.ArgumentParser(add_help=False)
        parser_vtam.add_argument('--db', **cls.args_db)
        parser_vtam.add_argument('--log', **cls.args_log_file)
        parser_vtam.add_argument('--threads', action='store',
                                     help="Number of threads",
                                     required=False,
                                     default=multiprocessing.cpu_count())
        parser_vtam.add_argument('-v', **cls.args_log_verbosity)
        parser_vtam.add_argument('--dry-run', '-n', dest='dryrun', action='store_true', required=False,
                                 help="Only display what would have been done.")
        parser_vtam.add_argument('-F', '--forceall', dest='forceall', action='store_true',
                                 help="Force argument of WopMars", required=False)
        parser_vtam.add_argument('--params', action='store', default=None, help="YML file with parameter values",
                                 required=False)
        parser_vtam.add_argument('-t', '--targetrule', dest='targetrule', action='store', default=None,
                                 help="Execute the workflow to the given target RULE: SampleInformation, ...",
                                 required=False)
        parser_vtam.add_argument('-f', '--sourcerule', dest='sourcerule', action='store', default=None,
                                 help="Execute the workflow from the given RULE.",
                                 required=False)
        subparsers = parser_vtam.add_subparsers()

        ################################################################################################################
        #
        # create the parser for the "merge" command
        #
        ################################################################################################################

        parser_vtam_merge = subparsers.add_parser('merge', add_help=True, parents=[parser_vtam])
        parser_vtam_merge.add_argument('--fastqinfo', action='store', help="TSV file with FASTQ sample information",
                                       required=True,
                                       type=lambda x: ArgParserChecker.check_file_exists_and_is_nonempty(x,
                                                                    error_message="Verify the '--fastqinfo' argument"))
        parser_vtam_merge\
            .add_argument('--fastainfo', action='store', help="REQUIRED: Output TSV file for FASTA sample information",
                          required=True)
        parser_vtam_merge.add_argument('--fastqdir', action='store', help="Directory with FASTQ files", required=True,
                                       type=lambda x:
                                       ArgParserChecker.check_dir_exists_and_is_nonempty(x,
                                                                    error_message="Verify the '--fastqdir' argument"))
        parser_vtam_merge.add_argument('--fastadir', action='store', help="Directory with FASTA files", required=True)
        parser_vtam_merge.set_defaults(command='merge')  # This attribute will trigget the good command

        parser_vtam_merge.add_argument('--outdir', action='store', help="REQUIRED: Directory for output", default="out",
                                     required=False)

        ################################################################################################################
        #
        # create the parser for the "asv" command
        #
        ################################################################################################################

        parser_vtam_asv = subparsers.add_parser('asv', add_help=True, parents=[parser_vtam])
        parser_vtam_asv\
            .add_argument('--fastainfo', action='store', help="REQUIRED: TSV file with FASTA sample information",
                          required=True, type=lambda x: ArgParserChecker
                          .check_file_exists_and_is_nonempty(x, error_message="Verify the '--fastainfo' argument"))
        parser_vtam_asv.add_argument('--fastadir', action='store', help="REQUIRED: Directory with FASTA files",
                                     required=True,
                                     type=lambda x:
                                     ArgParserChecker.check_file_exists_and_is_nonempty(x,
                                                                   error_message="Verify the '--fastadir' argument"))
        parser_vtam_asv.add_argument('--outdir', action='store', help="REQUIRED: Directory for output", default="out",
                                     required=True)

        parser_vtam_asv.add_argument('--blast_db', action='store',
                                     help="REQUIRED: Blast DB directory (Full or custom one) with nt files",
                                     required=True,
                                     type=lambda x: ArgParserChecker.check_blast_db_argument(x))
        parser_vtam_asv.add_argument('--taxonomy', dest='taxonomy', action='store',
                                     help="""REQUIRED: SQLITE DB with taxonomy information.

        This database is create with the command: vtam taxonomy. For instance

        vtam taxonomy -o taxonomy.sqlite to create a database in the current directory.""",
                                     required=True,
                                     type=lambda x: ArgParserChecker.check_file_exists_and_is_nonempty(x,
                                                                  error_message="Verify the '--taxonomy' argument"))
        parser_vtam_asv.add_argument('--threshold_specific', default=None, action='store', required=False,
                                 help="TSV file with variant (col1: variant; col2: threshold) or variant-replicate "
                                  "(col1: variant; col2: replicate; col3: threshold)specific thresholds. Header expected.",
                                 type=lambda x: ArgParserChecker.check_file_exists_and_is_nonempty(x,
                                              error_message="Verify the '--threshold_specific' argument"))
        parser_vtam_asv.set_defaults(command='asv')  # This attribute will trigget the good command

        ################################################################################################################
        #
        # create the parser for the "taxassign" command
        #
        ################################################################################################################

        parser_vtam_taxassign = subparsers.add_parser('taxassign', add_help=True)
        parser_vtam_taxassign.add_argument('--db', **cls.args_db)
        parser_vtam_taxassign.add_argument('--log', **cls.args_log_file)
        parser_vtam_taxassign.add_argument('-v', **cls.args_log_verbosity)
        parser_vtam_taxassign.add_argument('--blast_db', action='store',
                                     help="REQUIRED: Blast DB directory (Full or custom one) with nt files",
                                     required=True,
                                     type=lambda x: ArgParserChecker.check_blast_db_argument(x))
        parser_vtam_taxassign.add_argument('--taxonomy', dest='taxonomy', action='store',
                                     help="""REQUIRED: SQLITE DB with taxonomy information.

        This database is create with the command: vtam taxonomy. For instance

        vtam taxonomy -o taxonomy.sqlite to create a database in the current directory.""",
                                     required=True,
                                     type=lambda x: ArgParserChecker.check_file_exists_and_is_nonempty(x,
                                                                  error_message="Verify the '--taxonomy' argument"))
        parser_vtam_taxassign.set_defaults(command='taxassign')  # This attribute will trigget the good command

        ################################################################################################################
        #
        # create the parser for the "optimize" command
        #
        ################################################################################################################

        parser_vtam_optimize = subparsers.add_parser('optimize', add_help=True,  parents=[parser_vtam])
        parser_vtam_optimize\
            .add_argument('--fastainfo', action='store', help="REQUIRED: TSV file with FASTA sample information",
                          required=True, type=lambda x: ArgParserChecker
                          .check_file_exists_and_is_nonempty(x, error_message="Verify the '--fastainfo' argument"))
        parser_vtam_optimize.add_argument('--fastadir', action='store', help="REQUIRED: Directory with FASTA files",
                                     required=True,
                                     type=lambda x:
                                     ArgParserChecker.check_file_exists_and_is_nonempty(x,
                                                                   error_message="Verify the '--fastadir' argument"))
        parser_vtam_optimize.add_argument('--outdir', action='store', help="Directory for output", default="out",
                                          required=True)
        parser_vtam_optimize.add_argument('--variant_known', action='store', help="TSV file with known variants",
                                          required=True)
        parser_vtam_optimize.set_defaults(command='optimize')  # This attribute will trigget the good command

        ################################################################################################################
        #
        # create the parser for the "pool_markers" command
        #
        ################################################################################################################

        parser_vtam_pool_markers = subparsers.add_parser('pool_markers', add_help=True, formatter_class=argparse.RawTextHelpFormatter)
        parser_vtam_pool_markers.add_argument('--db', action='store', required=True, help="SQLITE file with DB")
        parser_vtam_pool_markers.add_argument('--runmarker', action='store', default=None,
                                     help="""Input TSV file with two columns and headers 'run_name' and 'marker_name'.
                                        Default: Uses all runs and markers in the DB
                                        Example:
                                        run_name	marker_name
                                        prerun	MFZR
                                        prerun	ZFZR""",
                                     required=False, type=ArgParserChecker.check_parser_vtam_pool_markers_arg_runmarker)
        parser_vtam_pool_markers.add_argument('--pooledmarkers', action='store', help="REQUIRED: Output TSV file with pooled markers",
                                       required=True)
        parser_vtam_pool_markers.set_defaults(command='pool_markers')  # This attribute will trigger the good command
        parser_vtam_pool_markers.add_argument('--taxonomy', dest='taxonomy', action='store',
                                     help="""REQUIRED: SQLITE DB with taxonomy information.

        This database is create with the command: vtam taxonomy. For instance

        vtam taxonomy -o taxonomy.sqlite to create a database in the current directory.""",
                                     required=True,
                                     type=lambda x: ArgParserChecker.check_file_exists_and_is_nonempty(x,
                                                                  error_message="Verify the '--taxonomy' argument"))

        ################################################################################################################
        #
        # create the parser for the "taxonomy" command
        #
        ################################################################################################################

        parser_vtam_taxonomy = subparsers.add_parser('taxonomy', add_help=True)
        parser_vtam_taxonomy.add_argument('-o', '--output', dest='output', action='store', help="Path to TSV taxonomy file",
                            required=True)
        parser_vtam_taxonomy.add_argument('--precomputed', dest='precomputed', action='store_true', default=False,
                            help="Will download precomputed taxonomy database, which is likely not the most recent one.",
                            required=False)
        parser_vtam_taxonomy.set_defaults(command='taxonomy')  # This attribute will trigger the good command

        ################################################################################################################
        #
        # create the parser for the "coi_db" command
        #
        ################################################################################################################

        parser_vtam_coi_blast_db = subparsers.add_parser('coi_blast_db', add_help=True)
        parser_vtam_coi_blast_db.add_argument('--coi_blast_db', dest='coi_blast_db', action='store', help="Path COI Blast DB",
                            required=True)
        parser_vtam_coi_blast_db.set_defaults(command='coi_blast_db')  # This attribute will trigger the good command

        return parser_vtam
