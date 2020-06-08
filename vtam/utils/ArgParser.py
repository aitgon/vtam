import argparse
import multiprocessing
import os
import pandas
import yaml

from vtam import CommandBlastCOI
from vtam.utils.CutoffSpecificFile import CutoffSpecificFile
from vtam.utils.KnownOccurrences import KnownOccurrences
from vtam.utils.ParamsFile import ParamsFile
from vtam.utils.SampleInformationFile import SampleInformationFile
from vtam.utils import constants
from vtam.utils.constants import header_merged_fasta, header_paired_fastq, header_sortedread_fasta, \
    coi_blast_db_gz_url


class ArgParserChecker(object):
    """Methods to check arguments"""

    @staticmethod
    def check_dir_exists_and_is_nonempty(path):
        """Checks if directory exists and is not empty

        :param path: Valid non-empty directory tsv_path
        :return: void
        """
        if not os.path.isdir(path):
            raise argparse.ArgumentTypeError(
                "The directory '{}' does not exist. Please fix it.".format(path))
        elif not len(os.listdir(path)) > 0:
            raise argparse.ArgumentTypeError(
                "The directory '{}' is empty. Please fix it.".format(path))
        else:
            return path

    @staticmethod
    def check_file_exists_and_is_nonempty(path):
        """Checks if file exists and is not empty

        :param path: Valid non-empty file tsv_path
        :return: void

        """
        if not os.path.isfile(path):
            raise argparse.ArgumentTypeError(
                "The file '{}' does not exist. Please fix it.".format(path))
        elif not os.stat(path).st_size > 0:
            raise argparse.ArgumentTypeError(
                "The file '{}' is empty. Please fix it.".format(path))
        else:
            return path  # return the tsv_path

    @classmethod
    def check_params_yml(cls, path):
        """Checks if file exists and is not empty

        :param path: Path to the sampleselect TSV file
        :return: void
        """

        path = cls.check_file_exists_and_is_nonempty(path)

        with open(path, 'r') as fin:
            params_user = yaml.load(fin, Loader=yaml.SafeLoader)
            for params_k in params_user:
                if not (params_k in constants.get_params_default_dic()):
                    raise argparse.ArgumentTypeError(
                        "This parameter '{}' in the YML file '{}' is not recognized by VTAM. "
                        "Please fix it.".format(
                            params_k, path))
                else:
                    return path

    @staticmethod
    def check_taxassign_taxonomy(path):
        """Check taxonomy format

        :param path: Valid non-empty file tsv_path
        :return: void

        """
        if not os.path.isfile(path):
            raise argparse.ArgumentTypeError(
                "The file '{}' does not exist. Please fix it.".format(path))
        elif not os.stat(path).st_size > 0:
            raise argparse.ArgumentTypeError(
                "The file '{}' is empty. Please fix it.".format(path))
        header_lower = {
            'tax_id',
            'parent_tax_id',
            'rank',
            'name_txt',
            'old_tax_id'}
        df = pandas.read_csv(path, sep="\t", header=0)
        df.columns = df.columns.str.lower()
        if set(
                df.columns) >= header_lower:  # contains at least the 'header_lower' columns
            return path
        else:
            raise argparse.ArgumentTypeError(
                "The format of file '{}' is wrong. Please fix it.".format(path))

    @staticmethod
    def check_taxassign_variants(path):
        """Check variants format

        :param path: Valid non-empty file tsv_path
        :return: void

        """
        if not os.path.isfile(path):
            raise argparse.ArgumentTypeError(
                "The file '{}' does not exist. Please fix it.".format(path))
        elif not os.stat(path).st_size > 0:
            raise argparse.ArgumentTypeError(
                "The file '{}' is empty. Please fix it.".format(path))
        header_lower = {'sequence'}
        df = pandas.read_csv(path, sep="\t", header=0)
        df.columns = df.columns.str.lower()
        if set(
                df.columns) >= header_lower:  # contains at least the 'header_lower' columns
            return path
        else:
            raise argparse.ArgumentTypeError(
                "The format of file '{}' is wrong. Please fix it.".format(path))


class ArgParser:

    args_db = {
        'dest': 'db',
        'action': 'store',
        'default': 'db.sqlite',
        'required': False,
        'help': "Database in SQLITE format"}
    args_log_file = {
        'dest': 'log_file',
        'action': 'store',
        'help': "Write log to file.",
        'required': False}
    args_log_verbosity = {
        'dest': 'log_verbosity',
        'action': 'count',
        'default': 0,
        'required': False,
        'help': "Set verbosity level, eg. None (Error level) -v (Info level)."}

    parser_vtam_main = None

    @classmethod
    def get_main_arg_parser(cls):
        """

        :return:
        """
        # create the top-level parser
        parser_vtam_main = argparse.ArgumentParser(add_help=False)
        parser_vtam_main.add_argument(
            '--params',
            action='store',
            default=None,
            help="YML file with parameter values",
            required=False,
            type=lambda x: ParamsFile(params_path=x).argparse_checker_params_file())
        parser_vtam_main.add_argument('--log', **cls.args_log_file)
        parser_vtam_main.add_argument('--threads', action='store',
                                      help="Number of threads",
                                      required=False,
                                      default=multiprocessing.cpu_count())
        parser_vtam_main.add_argument('-v', **cls.args_log_verbosity)
        subparsers = parser_vtam_main.add_subparsers()

        ############################################################################################
        #
        # create the parser for the "merge" command
        #
        ############################################################################################

        cls.create_merge(subparsers=subparsers, parent_parser=parser_vtam_main)

        #######################################################################
        #
        # create the parser for the "sortreads" command
        #
        #######################################################################

        cls.create_sortreads(subparsers=subparsers, parent_parser=parser_vtam_main)

        #######################################################################
        #
        # create the parser for the "filter" command
        #
        #######################################################################

        cls.create_filter(subparsers=subparsers, parent_parser=parser_vtam_main)

        #######################################################################
        #
        # create the parser for the "optimize" command
        #
        #######################################################################

        cls.create_optimize( subparsers=subparsers, parent_parser=parser_vtam_main)

        #######################################################################
        #
        # create the parser for the "pool" command
        #
        #######################################################################

        cls.create_pool(subparsers=subparsers, parent_parser=parser_vtam_main)

        #######################################################################
        #
        # create the parser for the "taxassign" command
        #
        #######################################################################

        cls.create_taxassign(subparsers=subparsers, parent_parser=parser_vtam_main)

        #######################################################################
        #
        # create the parser for the "taxonomy" command
        #
        #######################################################################

        cls.create_taxonomy(subparsers=subparsers, parent_parser=parser_vtam_main)

        #######################################################################
        #
        # create the parser for the "coi_db" command
        #
        #######################################################################

        cls.create_coiblastdb(subparsers=subparsers)

        return parser_vtam_main

    @classmethod
    def create_merge(cls, subparsers, parent_parser):

        parser_vtam_merge = subparsers.add_parser('merge', add_help=True,
            formatter_class=argparse.RawTextHelpFormatter,
            parents=[parent_parser])

        parser_vtam_merge.add_argument('--fastqinfo', action='store',
            help="TSV file with FASTQ sample information",
            required=True,
            type=lambda x: SampleInformationFile(x).check_args(
                header=header_paired_fastq))

        parser_vtam_merge .add_argument( '--fastainfo',
            action='store',
            help="REQUIRED: Output TSV file for FASTA sample information",
            required=True)

        parser_vtam_merge.add_argument(
            '--fastqdir',
            action='store',
            help="Directory with FASTQ files",
            required=True,
            type=ArgParserChecker.check_dir_exists_and_is_nonempty)

        parser_vtam_merge.add_argument(
            '--fastadir',
            action='store',
            help="Directory with FASTA files",
            required=True)
        # This attribute will trigger the good command

        parser_vtam_merge.set_defaults(command='merge')

    @classmethod
    def create_sortreads(cls, subparsers, parent_parser):

        parser_vtam_sortreads = subparsers.add_parser(
            'sortreads',
            add_help=True,
            formatter_class=argparse.RawTextHelpFormatter,
            parents=[parent_parser])

        parser_vtam_sortreads .add_argument(
            '--fastainfo',
            action='store',
            help="REQUIRED: TSV file with FASTA information",
            required=True,
            type=lambda x: SampleInformationFile(x).check_args(
                header=header_merged_fasta))

        parser_vtam_sortreads.add_argument(
            '--fastadir',
            action='store',
            help="REQUIRED: Directory with FASTA files",
            required=True,
            type=ArgParserChecker.check_dir_exists_and_is_nonempty)

        parser_vtam_sortreads.add_argument(
            '--outdir',
            action='store',
            help="REQUIRED: Output directory for trimmed and demultiplexed files",
            default="out",
            required=True)
        # This attribute will trigger the good command

        parser_vtam_sortreads.set_defaults(command='sortreads')

    @classmethod
    def create_filter(cls, subparsers, parent_parser):

        parser_vtam_filter = subparsers.add_parser(
            'filter', add_help=True, parents=[parent_parser])

        parser_vtam_filter .add_argument(
            '--readinfo',
            action='store',
            help="REQUIRED: TSV file with information of sorted read files",
            required=True,
            type=lambda x: SampleInformationFile(x).check_args(
                header=header_sortedread_fasta)
        )
        parser_vtam_filter.add_argument(
            '--readdir',
            action='store',
            help="REQUIRED: TSV file with information of sorted read files",
            required=True,
            type=ArgParserChecker.check_dir_exists_and_is_nonempty
        )
        parser_vtam_filter .add_argument(
            '--asvtable',
            action='store',
            help="REQUIRED: Output TSV file for the amplicon sequence variants (ASV) table",
            required=True)

        parser_vtam_filter.add_argument(
            '--cutoff_specific',
            dest='cutoff_specific',
            default=None,
            action='store',
            required=False,
            help="TSV file with variant (col1: variant; col2: cutoff) or variant-replicate "
            "(col1: variant; col2: replicate; col3: cutoff)specific cutoffs.",
            type=lambda x: CutoffSpecificFile(x).argparse_checker())

        parser_vtam_filter.add_argument(
            '--lfn_variant_replicate',
            action='store_true',
            help="If set, VTAM will run the algorithm for the low frequency noise over variant and replicates",
            required=False,
            default=False)

        #######################################################################
        #
        # Wopmars args
        #
        #######################################################################

        parser_vtam_filter.add_argument('--db', **cls.args_db)

        parser_vtam_filter.add_argument(
            '--dry-run_name',
            '-n',
            dest='dryrun',
            action='store_true',
            required=False,
            help="Only display what would have been done.")

        parser_vtam_filter.add_argument(
            '-F',
            '--forceall',
            dest='forceall',
            action='store_true',
            help="Force argument of WopMars",
            required=False)

        parser_vtam_filter.add_argument(
            '-U',
            '--until',
            dest='until',
            action='store',
            default=None,
            help="Execute the workflow until the given target RULE: SampleInformation, ...",
            required=False)

        parser_vtam_filter.add_argument(
            '-S',
            '--since',
            dest='since',
            action='store',
            default=None,
            help="Execute the workflow since the given RULE.",
            required=False)

        # This attribute will trigger the good command
        parser_vtam_filter.set_defaults(command='filter')

    @classmethod
    def create_optimize(cls, subparsers, parent_parser):

        parser_vtam_optimize = subparsers.add_parser(
            'optimize', add_help=True, parents=[parent_parser])

        parser_vtam_optimize .add_argument(
            '--readinfo',
            action='store',
            help="REQUIRED: TSV file with information of sorted read files",
            required=True,
            type=lambda x: SampleInformationFile(x).check_args(
                header=header_sortedread_fasta))

        parser_vtam_optimize.add_argument(
            '--readdir',
            action='store',
            help="REQUIRED: Directory with sorted read files",
            required=True,
            type=ArgParserChecker.check_dir_exists_and_is_nonempty)

        parser_vtam_optimize.add_argument(
            '--outdir',
            action='store',
            help="Directory for output",
            default="out",
            required=True)

        parser_vtam_optimize.add_argument(
            '--known_occurrences',
            action='store',
            help="TSV file with known variants",
            required=True,
            type=lambda x: KnownOccurrences(x).argparse_checker_known_occurrences())

        parser_vtam_optimize.add_argument(
            '--lfn_variant_replicate',
            action='store_true',
            help="If set, VTAM will run the algorithm for the low frequency noise over variant and replicates",
            required=False,
            default=False)

        #######################################################################
        #
        # Wopmars args
        #
        #######################################################################

        parser_vtam_optimize.add_argument('--db', **cls.args_db)

        parser_vtam_optimize.add_argument(
            '--dry-run_name',
            '-n',
            dest='dryrun',
            action='store_true',
            required=False,
            help="Only display what would have been done.")

        parser_vtam_optimize.add_argument(
            '-F',
            '--forceall',
            dest='forceall',
            action='store_true',
            help="Force argument of WopMars",
            required=False)

        parser_vtam_optimize.add_argument(
            '-U',
            '--until',
            dest='until',
            action='store',
            default=None,
            help="Execute the workflow until the given target RULE: SampleInformation, ...",
            required=False)
        parser_vtam_optimize.add_argument(
            '-S',
            '--since',
            dest='since',
            action='store',
            default=None,
            help="Execute the workflow since the given RULE.",
            required=False)

        # This attribute will trigger the good command
        parser_vtam_optimize.set_defaults(command='optimize')

    @classmethod
    def create_pool(cls, subparsers, parent_parser):

        parser_vtam_pool_markers = subparsers.add_parser(
            'pool',
            add_help=True,
            formatter_class=argparse.RawTextHelpFormatter,
            parents=[parent_parser])
        parser_vtam_pool_markers.add_argument(
            '--db', action='store', required=True, help="SQLITE file with DB")
        from vtam.utils.RunMarkerFile import RunMarkerFile
        parser_vtam_pool_markers.add_argument(
            '--runmarker',
            action='store',
            default=None,
            help=RunMarkerFile.help(),
            required=True,
            type=lambda x: RunMarkerFile(x).check_argument())
        parser_vtam_pool_markers.add_argument(
            '--output',
            action='store',
            help="REQUIRED: Output TSV file with pooled markers",
            required=True)

        # This attribute will trigger the good command
        parser_vtam_pool_markers.set_defaults(command='pool')

    @classmethod
    def create_taxassign(cls, subparsers, parent_parser):

        parser_vtam_taxassign = subparsers.add_parser(
            'taxassign',
            add_help=True,
            formatter_class=argparse.RawTextHelpFormatter,
            parents=[parent_parser])

        parser_vtam_taxassign .add_argument(
            '--variants',
            action='store',
            help="REQUIRED: TSV file with variant sequences and sequence header in the last column.",
            required=True,
            type=lambda x: ArgParserChecker.check_taxassign_variants(x))
        parser_vtam_taxassign .add_argument(
            '--output',
            action='store',
            help="REQUIRED: TSV file where the taxon assignation has beeen added.",
            required=True)
        parser_vtam_taxassign.add_argument(
            '--mode',
            dest='mode',
            default="unassigned",
            action='store',
            required=False,
            choices=[
                'unassigned',
                'reset'],
            help="The default 'unassigned' mode will only assign 'unassigned' variants."
            "The alternative 'reset' mode will erase the TaxAssign table and reassigned all "
            "input variants.")
        parser_vtam_taxassign.add_argument('--db', **cls.args_db)
        parser_vtam_taxassign.add_argument(
            '--blastdbdir',
            action='store',
            help="REQUIRED: Blast DB directory (Full or custom one) with DB files.",
            required=True,
            type=ArgParserChecker.check_dir_exists_and_is_nonempty)
        parser_vtam_taxassign.add_argument(
            '--blastdbname',
            action='store',
            help="REQUIRED: Blast DB name. It corresponds to file name (without suffix)"
            "of blast DB files.",
            required=True)
        parser_vtam_taxassign.add_argument(
            '--taxonomy',
            dest='taxonomy',
            action='store',
            help="""REQUIRED: SQLITE DB with taxonomy information.

        This database is create with the command: vtam taxonomy. For instance

        vtam taxonomy -o taxonomy.sqlite to create a database in the current directory.""",
            required=True,
            type=ArgParserChecker.check_taxassign_taxonomy)

        # This attribute will trigger the good command
        parser_vtam_taxassign.set_defaults(command='taxassign')

    @classmethod
    def create_taxonomy(cls, subparsers, parent_parser):

        parser_vtam_taxonomy = subparsers.add_parser('taxonomy', add_help=True,
                                                     parents=[parent_parser])
        parser_vtam_taxonomy.add_argument(
            '-o',
            '--output',
            dest='output',
            action='store',
            help="Path to TSV taxonomy file",
            required=True)
        parser_vtam_taxonomy.add_argument(
            '--precomputed',
            dest='precomputed',
            action='store_true',
            default=False,
            help="Will download precomputed taxonomy database, "
                 "which is likely not the most recent one.",
            required=False)
        # This attribute will trigger the good command
        parser_vtam_taxonomy.set_defaults(command='taxonomy')

    @classmethod
    def create_coiblastdb(cls, subparsers):

        parser_vtam_coi_blast_db = subparsers.add_parser(
            'coi_blast_db', add_help=True)
        parser_vtam_coi_blast_db.add_argument(
            '--blastdbdir',
            dest='blastdbdir',
            action='store',
            help="Path COI Blast DB",
            required=True)
        parser_vtam_coi_blast_db.add_argument(
            '--blastdbname',
            dest='blastdbname',
            action='store',
            help="COI Blast DB name, eg coi_blast_db_20191211, coi_blast_db_20200420. Versions can be found as basename here: {}".format(os.path.dirname(coi_blast_db_gz_url)),
            required=False,
            default='coi_blast_db',
            type=lambda x: CommandBlastCOI(x).argparse_checker_blast_coi_blastdbname(),
        )
        # This attribute will trigger the good command
        parser_vtam_coi_blast_db.set_defaults(command='coi_blast_db')
