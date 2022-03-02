import argparse
import multiprocessing
import os
import pathlib
import pandas
import yaml

from vtam import CommandBlastCOI
from vtam.utils.FileCutoffSpecific import FileCutoffSpecific
from vtam.utils.FileKnownOccurrences import FileKnownOccurrences
from vtam.utils.FileParams import FileParams
from vtam.utils.FileSampleInformation import FileSampleInformation
from vtam.utils import constants
from vtam.utils.constants import header_merged_fasta, header_paired_fastq, header_sortedread_fasta, \
    coi_blast_db_gz_url1
import vtam


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
                "The --asvtable TSV file must contain a column with a 'sequence' header: '{}'.".format(
                    path))


class ArgParser:

    ############################################################################################
    #
    # Specific parsers
    #
    ############################################################################################

    parser_params = argparse.ArgumentParser(add_help=False)
    parser_params.add_argument(
        '--params', action='store', default=None, help="YML file with parameter values",
        required=False, type=lambda x: FileParams(
            params_path=x).argparse_checker_params_file())

    parser_log = argparse.ArgumentParser(add_help=False)
    parser_log.add_argument('--log', dest='log', action='store', help="write log to LOG file.",
                            required=False)

    parser_threads = argparse.ArgumentParser(add_help=False)
    parser_threads.add_argument('--threads', dest='threads', action='store', help="number of threads",
                                  required=False, default=multiprocessing.cpu_count())

    parser_verbosity = argparse.ArgumentParser(add_help=False)
    parser_verbosity.add_argument('-v', dest='log_verbosity', action='count', default=0,
                                  required=False, help="set verbosity level -v or -vv")

    parser_wopmars_db = argparse.ArgumentParser(add_help=False)
    parser_wopmars_db.add_argument('--db', dest='db', action='store', default='db.sqlite',
                                   required=False, help="database file in SQLITE format")

    parser_wopmars_dryrun = argparse.ArgumentParser(add_help=False)
    parser_wopmars_dryrun.add_argument(
        '--dry-run',
        '-n',
        dest='dryrun',
        action='store_true',
        required=False,
        help="displays only command out without running it")

    parser_wopmars_forceall = argparse.ArgumentParser(add_help=False)
    parser_wopmars_forceall.add_argument(
        '-F',
        '--forceall',
        dest='forceall',
        action='store_true',
        help="force rerun all rules",
        required=False)

    parser_vtam_main = None

    @classmethod
    def get_main_arg_parser(cls):
        """

        :return:
        """

        ############################################################################################
        #
        # Top-level parser
        #
        ############################################################################################

        # config = RawConfigParser()
        # config.read(os.path.join(PathManager.get_package_path(), 'setup.cfg'))
        # version = config.get('metadata', 'version')

        parser_vtam_main = argparse.ArgumentParser(prog='vtam', description='%(prog)s {} - VTAM - Validation and Taxonomic Assignation of Metabarcoding Data'.format(vtam.__version__))
        parser_vtam_main.add_argument('--version', action='version', version='%(prog)s {}'.format(vtam.__version__))
        subparsers = parser_vtam_main.add_subparsers(title='VTAM sub-commands')

        ############################################################################################
        #
        # create the parsers
        #
        ############################################################################################

        cls.add_parser_example(subparsers=subparsers)

        cls.add_parser_merge(subparsers=subparsers)

        cls.add_parser_sortreads(subparsers=subparsers)

        cls.add_parser_filter(subparsers=subparsers)

        cls.add_parser_optimize(subparsers=subparsers)

        cls.add_parser_pool(subparsers=subparsers)

        cls.add_parser_taxassign(subparsers=subparsers)

        cls.add_parser_taxonomy(subparsers=subparsers)

        cls.add_parser_coiblastdb(subparsers=subparsers)

        return parser_vtam_main

    @classmethod
    def add_parser_example(cls, subparsers):
        parser_vtam_merge = subparsers.add_parser('example', add_help=True,
                                                  parents=[cls.parser_params, cls.parser_log,
                                                           cls.parser_threads, cls.parser_verbosity],
                                                  help="generates data for quick start")

        parser_vtam_merge.add_argument('--outdir', action='store',
                                       help="directory for quick start data",
                                       required=False, default='example', type=lambda x: pathlib.Path(x).mkdir(exist_ok=True, parents=True) or x)

        parser_vtam_merge.set_defaults(command='example')

    @classmethod
    def add_parser_merge(cls, subparsers):
        parser_vtam_merge = subparsers.add_parser('merge', add_help=True,
                                                  parents=[cls.parser_params, cls.parser_log,
                                                           cls.parser_threads, cls.parser_verbosity],
                                                  help="merges paired-end reads")

        parser_vtam_merge.add_argument('--fastqinfo', action='store',
                                       help="input TSV file with paired FASTQ file information",
                                       required=True,
                                       type=lambda x: FileSampleInformation(x).check_args(
                                           header=header_paired_fastq))

        parser_vtam_merge.add_argument('--fastainfo',
                                       action='store',
                                       help="output TSV file with merged FASTA file information",
                                       required=True)

        parser_vtam_merge.add_argument(
            '--fastqdir',
            action='store',
            help="input directory with paired FASTQ files",
            required=True,
            type=ArgParserChecker.check_dir_exists_and_is_nonempty)

        parser_vtam_merge.add_argument(
            '--fastadir',
            action='store',
            help="output directory with merged FASTA files",
            required=True)
        # This attribute will trigger the good command

        parser_vtam_merge.set_defaults(command='merge')

    @classmethod
    def add_parser_sortreads(cls, subparsers):
        parser_vtam_sortreads = subparsers.add_parser(
            'sortreads',
            add_help=True,
            parents=[cls.parser_params, cls.parser_log,
            cls.parser_threads, cls.parser_verbosity],
            help="sorts (Trims and demultiplexes) reads to biological samples and replicates according to the presence of sequence tags and primers")

        parser_vtam_sortreads.add_argument(
            '--fastainfo',
            action='store',
            help="input TSV file with FASTA file information",
            required=True,
            type=lambda x: FileSampleInformation(x).check_args(
                header=header_merged_fasta))

        parser_vtam_sortreads.add_argument(
            '--fastadir',
            action='store',
            help="input directory with FASTA files",
            required=True,
            type=ArgParserChecker.check_dir_exists_and_is_nonempty)

        parser_vtam_sortreads.add_argument(
            '--sorteddir',
            action='store',
            help="output directory with sorted reads (Trimmed and demultiplexed) in FASTA files and TSV file with corresponnding FASTA file information ('SORTEDDIR/sortedinfo.tsv')",
            default="out",
            required=True)
        # This attribute will trigger the good command

        parser_vtam_sortreads.set_defaults(command='sortreads')

    @classmethod
    def add_parser_filter(cls, subparsers):
        parser_vtam_filter = subparsers.add_parser(
            'filter', add_help=True,
            parents=[cls.parser_params, cls.parser_log,
            cls.parser_threads, cls.parser_verbosity, cls.parser_wopmars_db,
                     cls.parser_wopmars_dryrun, cls.parser_wopmars_forceall],
            help="filters out sequence artifacts and creates an amplicon sequence variant (ASV) table.")

        parser_vtam_filter.add_argument(
            '--sortedinfo',
            action='store',
            help="input TSV file with information about FASTA files containing sorted reads",
            required=True,
            type=lambda x: FileSampleInformation(x).check_args(
                header=header_sortedread_fasta)
        )
        parser_vtam_filter.add_argument(
            '--sorteddir',
            action='store',
            help="input directory with sorted (Trimmed and demultiplexed) FASTA files",
            required=True,
            type=ArgParserChecker.check_dir_exists_and_is_nonempty
        )
        parser_vtam_filter.add_argument(
            '--asvtable',
            action='store',
            help="output TSV file for the amplicon sequence variants (ASV) table",
            required=True)

        parser_vtam_filter.add_argument('--cutoff_specific', dest='cutoff_specific', default=None,
                                        action='store',
                                        required=False,
                                        help="TSV file with variant (col1: variant; col2: cutoff) or variant-replicate "
                                             "(col1: variant; col2: replicate; col3: cutoff)specific cutoffs",
                                        type=lambda x: FileCutoffSpecific(x).argparse_checker())

        parser_vtam_filter.add_argument(
            '--lfn_variant_replicate',
            action='store_true',
            help="if set, VTAM will run the algorithm for the low frequency noise over variant and replicates",
            required=False,
            default=False)

        parser_vtam_filter.add_argument(
            '--known_occurrences',
            action='store',
            help="TSV file with expected (keep) occurrences",
            required=False,
            type=lambda x: FileKnownOccurrences(x).argparse_checker_known_occurrences())

        parser_vtam_filter.add_argument(
            '-U',
            '--until',
            dest='until',
            action='store',
            default=None,
            help="""execute '%(prog)s' UNTIL one rule, where the rule order looks like:            
1. SampleInformation, 2. VariantReadCount, 3. FilterLFN, 4. FilterMinReplicateNumber, 5. FilterPCRerror, 6. FilterChimera, 7. FilterMinReplicateNumber2, 8. FilterRenkonen, 9. FilterMinReplicateNumber3, 10. FilterIndel, 11. FilterCodonStop, 12. ReadCountAverageOverReplicates, 13. MakeAsvTable""",
            required=False)

        parser_vtam_filter.add_argument(
            '-S',
            '--since',
            dest='since',
            action='store',
            default=None,
            help="""execute '%(prog)s' SINCE one rule, where the rule order looks like:
            1. SampleInformation, 2. VariantReadCount, 3. FilterLFN, 4. FilterMinReplicateNumber, 5. FilterPCRerror, 6. FilterChimera, 7. FilterMinReplicateNumber2, 8. FilterRenkonen, 9. FilterMinReplicateNumber3, 10. FilterIndel, 11. FilterCodonStop, 12. ReadCountAverageOverReplicates, 13. MakeAsvTable""",
            required=False)

        # This attribute will trigger the good command
        parser_vtam_filter.set_defaults(command='filter')

    @classmethod
    def add_parser_taxassign(cls, subparsers):
        parser_vtam_taxassign = subparsers.add_parser(
            'taxassign',
            add_help=True, parents=[cls.parser_params, cls.parser_log,
                                                           cls.parser_threads, cls.parser_verbosity,
                                    cls.parser_wopmars_db],
            help="assigns amplicon sequence variants (ASVs) to taxonomic groups")

        parser_vtam_taxassign.add_argument(
            '--asvtable',
            action='store',
            help="input TSV file with variant sequences and sequence header in the last column",
            required=True,
            type=lambda x: ArgParserChecker.check_taxassign_variants(x))
        parser_vtam_taxassign.add_argument(
            '--output',
            action='store',
            help="output TSV file where the assigned taxa have been added",
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
            help="the default 'unassigned' mode will only assign 'unassigned' variants"
                 "The alternative 'reset' mode will erase the TaxAssign table and reassigned all "
                 "input variants")
        parser_vtam_taxassign.add_argument(
            '--blastdbdir',
            action='store',
            help="input directory with (Full or custom one) Blast database files",
            required=True,
            type=ArgParserChecker.check_dir_exists_and_is_nonempty)
        parser_vtam_taxassign.add_argument(
            '--blastdbname',
            action='store',
            help="input Blast database name, which corresponds to the file name without suffix of the Blast database files",
            required=True)
        parser_vtam_taxassign.add_argument(
            '--taxonomy',
            dest='taxonomy',
            action='store',
            help="""input TSV file with taxonomy information.
        This file is created with the 'taxonomy' sub-command. For instance
        'vtam taxonomy -o taxonomy.tsv' creates the 'taxonomy.tsv' file in the current directory""",
            required=True,
            type=ArgParserChecker.check_taxassign_taxonomy)

        # This attribute will trigger the good command
        parser_vtam_taxassign.set_defaults(command='taxassign')

    @classmethod
    def add_parser_optimize(cls, subparsers):
        parser_vtam_optimize = subparsers.add_parser(
            'optimize', add_help=True,
            parents=[cls.parser_params, cls.parser_log,
            cls.parser_threads, cls.parser_verbosity, cls.parser_wopmars_db,
                     cls.parser_wopmars_dryrun, cls.parser_wopmars_forceall],
            help="finds out optimal parameters for filtering")

        parser_vtam_optimize.add_argument(
            '--sortedinfo',
            action='store',
            help="input TSV file with information about FASTA files containing sorted (trimmed and demultiplexed) reads",
            required=True,
            type=lambda x: FileSampleInformation(x).check_args(
                header=header_sortedread_fasta))

        parser_vtam_optimize.add_argument(
            '--sorteddir',
            action='store',
            help="input directory with sorted (Trimmed and demultiplexed) FASTA files",
            required=True,
            type=ArgParserChecker.check_dir_exists_and_is_nonempty)

        parser_vtam_optimize.add_argument(
            '-o',
            '--outdir',
            action='store',
            help="output directory",
            default="out",
            required=True)

        parser_vtam_optimize.add_argument(
            '--known_occurrences',
            action='store',
            help="TSV file with known variants",
            required=True,
            type=lambda x: FileKnownOccurrences(x).argparse_checker_known_occurrences())

        parser_vtam_optimize.add_argument(
            '--lfn_variant_replicate',
            action='store_true',
            help="if set, VTAM will run the algorithm for the low frequency noise over variant and replicates",
            required=False,
            default=False)

        parser_vtam_optimize.add_argument(
            '-U',
            '--until',
            dest='until',
            action='store',
            default=None,
            help="""executes '%(prog)s' UNTIL one rule, where the rules follow this order:
            1. SampleInformation, 2. VariantReadCount, 3. either OptimizeLFNsampleReplicate or OptimizePCRerror or OptimizeLFNreadCountAndLFNvariant""",
            required=False)
        parser_vtam_optimize.add_argument(
            '-S',
            '--since',
            dest='since',
            action='store',
            default=None,
            help="""executes '%(prog)s' SINCE one rule, where the rules follow this order: 
            1. SampleInformation, 2. VariantReadCount, 3. either OptimizeLFNsampleReplicate or OptimizePCRerror or OptimizeLFNreadCountAndLFNvariant""",
            required=False)

        # This attribute will trigger the good command
        parser_vtam_optimize.set_defaults(command='optimize')

    @classmethod
    def add_parser_pool(cls, subparsers):
        parser_vtam_pool_markers = subparsers.add_parser(
            'pool',
            add_help=True,
                                                  parents=[cls.parser_params, cls.parser_log,
                                                           cls.parser_threads, cls.parser_verbosity],
            help="pools amplicon sequence variants (ASVs) from different but overlapping markers")

        parser_vtam_pool_markers.add_argument(
            '--db', action='store', required=True, help="SQLITE file with DB")

        from vtam.utils.FileRunMarker import FileRunMarker
        parser_vtam_pool_markers.add_argument(
            '--runmarker',
            action='store',
            default=None,
            help=FileRunMarker.help(),
            required=True,
            type=lambda x: FileRunMarker(x).check_argument())

        parser_vtam_pool_markers.add_argument(
            '--asvtable',
            action='store',
            help="output TSV file with pooled markers and their occurrences in biological samples",
            required=True)

        parser_vtam_pool_markers.add_argument(
            '--readcounts',
            action='store_true',
            help="Default: False. If False, presence/absence of reads in sample is given."
                 "If True, sum of reads over pooled runs et/ou markers is given",
            required=False,
            default=False)

        # This attribute will trigger the good command
        parser_vtam_pool_markers.set_defaults(command='pool')

    @classmethod
    def add_parser_taxonomy(cls, subparsers):
        parser_vtam_taxonomy = subparsers.add_parser('taxonomy', add_help=True,
                                                  parents=[],
                                                     help="downloads a TSV file with the NCBI taxonomy information")

        parser_vtam_taxonomy.add_argument(
            '-o',
            '--output',
            dest='output',
            action='store',
            help="default: taxonomy.tsv. Path to TSV taxonomy file",
            required=False,
        default=os.path.join(os.getcwd(), 'taxonomy.tsv'))
        parser_vtam_taxonomy.add_argument(
            '--precomputed',
            dest='precomputed',
            action='store_true',
            default=False,
            help="default: False. Downloads precomputed taxonomy database, "
                 "which is likely an older database",
            required=False)
        # This attribute will trigger the good command
        parser_vtam_taxonomy.set_defaults(command='taxonomy')

    @classmethod
    def add_parser_coiblastdb(cls, subparsers):
        parser_vtam_coi_blast_db = subparsers.add_parser(
            'coi_blast_db', add_help=True,
            help="downloads a precomputed BLAST database for the cytochrome C oxidase subunit I (COI) marker")

        parser_vtam_coi_blast_db.add_argument(
            '--blastdbdir',
            dest='blastdbdir',
            action='store',
            help="output directory with custom Blast database files of the cytochrome C oxidase subunit I (COI) marker files",
            required=False,
            default='blastdb')
        parser_vtam_coi_blast_db.add_argument(
            '--blastdbname',
            dest='blastdbname',
            action='store',
            help="cytochrome C oxidase subunit I (COI) Blast database name among these current possibilities: coi_blast_db, coi_blast_db_20200420. Other versions if available can be found here: {}".format(
                os.path.dirname(coi_blast_db_gz_url1)),
            required=False,
            default='coi_blast_db',
            type=lambda x: CommandBlastCOI(x).argparse_checker_blast_coi_blastdbname(),
        )
        # This attribute will trigger the good command
        parser_vtam_coi_blast_db.set_defaults(command='coi_blast_db')
