import os
import sys

from vtam.CommandFilterOptimize import CommandFilterOptimize
from vtam.CommandBlastCOI import CommandBlastCOI
from vtam.CommandMerge import CommandMerge
from vtam.CommandPoolRunMarkers import CommandPoolRunMarkers
from vtam.CommandSortReads import CommandSortReads
from vtam.CommandTaxAssign import CommandTaxAssign
from vtam.CommandTaxonomy import CommandTaxonomy
from vtam.utils.ArgParser import ArgParser
from vtam.utils.Logger import Logger
from vtam.utils.Logger import LoggerArguments
from vtam.utils.VTAMexception import VTAMexception
from vtam.utils.WopmarsRunner import WopmarsRunner
from vtam.utils.constants import FilterLFNreference_records


class VTAM(object):

    usage_message = """usage: vtam <command> [<args>]

        These are the VTAM commands:

   merge      Merges paired-end reads
   sortreads  Trims, demultiplex and and sorts reads
   filter     Filters out sequence artifacts and creates an amplicon sequence variant (ASV) table
   taxassign     Assigns amplicon sequence variants to taxa
   optimize   Finds optimal parameters for filtering
   pool   Pools overlapping markers from the ASV table into one
   taxonomy   Creates the taxonomy TSV file required for tax assignation
   coi_blast_db   Downloads a precomputed COI Blast database
"""

    def __init__(self, sys_argv):

        #######################################################################
        #
        # Parse arguments
        #
        #######################################################################

        self.sys_argv = sys_argv
        # AG do not use abspath for the moment. Maybe later it can be used as
        # option
        parser = ArgParser.get_main_arg_parser()
        self.args = parser.parse_args(sys_argv)

        arg_parser_dic = vars(self.args)

        #######################################################################
        #
        # Parse log arguments
        #
        #######################################################################

        (LoggerArguments.instance()).update({'log_verbosity': arg_parser_dic['log_verbosity'],
                                             'log_file': arg_parser_dic['log_file']})
        if 'log_verbosity' in arg_parser_dic:
            os.environ['VTAM_LOG_VERBOSITY'] = str(
                arg_parser_dic['log_verbosity'])
        if 'log_file' in arg_parser_dic:
            os.environ['VTAM_LOG_FILE'] = str(arg_parser_dic['log_file'])

        #######################################################################
        #
        # Set arguments, logger
        #
        #######################################################################

        # Some arguments will be passed through environmental variables
        if 'threads' in arg_parser_dic:
            os.environ['VTAM_THREADS'] = str(arg_parser_dic['threads'])

        ###############################################################
        #
        # Subcommands: wopfile-dependent, filter, optimize
        #
        ###############################################################

        if arg_parser_dic['command'] in ['filter', 'optimize']:

            CommandFilterOptimize.main(arg_parser_dic=arg_parser_dic)

        ############################################################################################
        #
        # Subcommand: merge
        #
        ############################################################################################

        elif arg_parser_dic['command'] == 'merge':
            fastqinfo = arg_parser_dic['fastqinfo']
            fastqdir = arg_parser_dic['fastqdir']
            fastainfo = arg_parser_dic['fastainfo']
            fastadir = arg_parser_dic['fastadir']
            num_threads = arg_parser_dic['threads']
            params = arg_parser_dic['params']
            CommandMerge.main(fastqinfo=fastqinfo, fastqdir=fastqdir, fastainfo=fastainfo,
                              fastadir=fastadir, params=params, num_threads=num_threads)

        ############################################################################################
        #
        # Subcommand: sortreads
        #
        ############################################################################################

        elif arg_parser_dic['command'] == 'sortreads':
            fastadir = arg_parser_dic['fastadir']
            fastainfo = arg_parser_dic['fastainfo']
            outdir = arg_parser_dic['outdir']
            num_threads = arg_parser_dic['threads']
            params = arg_parser_dic['params']
            CommandSortReads.main(fastainfo=fastainfo, fastadir=fastadir, params=params,
                                  num_threads=num_threads, outdir=outdir)

        ############################################################################################
        #
        # Subcommand: filter
        #
        ############################################################################################

        elif arg_parser_dic['command'] == 'sortreads':
            fastadir = arg_parser_dic['fastadir']
            fastainfo = arg_parser_dic['fastainfo']
            outdir = arg_parser_dic['outdir']
            num_threads = arg_parser_dic['threads']
            params = arg_parser_dic['params']
            CommandSortReads.main(fastainfo=fastainfo, fastadir=fastadir, params=params,
                                  num_threads=num_threads, outdir=outdir)

        ############################################################################################
        #
        # Subcommand: taxassign
        #
        ############################################################################################

        elif arg_parser_dic['command'] == 'taxassign':
            db = arg_parser_dic['db']
            variants_tsv = arg_parser_dic['variants']
            output = arg_parser_dic['output']
            mode = arg_parser_dic['mode']
            taxonomy_tsv = arg_parser_dic['taxonomy']
            blasdb_dir_path = arg_parser_dic['blastdbdir']
            blastdbname_str = arg_parser_dic['blastdbname']
            num_threads = arg_parser_dic['threads']
            params = arg_parser_dic['params']
            CommandTaxAssign.main(db=db, mode=mode, variants_tsv=variants_tsv, output=output,
                taxonomy_tsv=taxonomy_tsv, blastdb_dir_path=blasdb_dir_path,
                blastdbname_str=blastdbname_str, params=params, num_threads=num_threads)

        ############################################################################################
        #
        # Subcommand: pool
        #
        ############################################################################################

        elif arg_parser_dic['command'] == 'pool':
            db = arg_parser_dic['db']
            run_marker_tsv = arg_parser_dic['runmarker']
            pooled_marker_tsv = arg_parser_dic['output']
            CommandPoolRunMarkers.main(db=db, pooled_marker_tsv=pooled_marker_tsv,
                run_marker_tsv=run_marker_tsv)

        ############################################################################################
        #
        # Subcommand: taxonomy
        #
        ############################################################################################

        elif arg_parser_dic['command'] == 'taxonomy':
            taxonomy_tsv = arg_parser_dic['output']
            precomputed = arg_parser_dic['precomputed']
            taxonomy = CommandTaxonomy(taxonomy_tsv=taxonomy_tsv)
            taxonomy.main(precomputed=precomputed)

        ############################################################################################
        #
        # Subcommand: coi blast
        #
        ############################################################################################

        elif arg_parser_dic['command'] == 'coi_blast_db':
            coi_blast_db_dir = arg_parser_dic['coi_blast_db_dir']
            coi_blast_db = CommandBlastCOI(coi_blast_db_dir=coi_blast_db_dir)
            coi_blast_db.download()

        ############################################################################################
        #
        # Else: run_name usage message
        #
        ############################################################################################

        else:
            print(VTAM.usage_message)


def main():
    if sys.argv[1:] == []:  # No arguments
        print(VTAM.usage_message)
        sys.exit(1)
    VTAM(sys.argv[1:])
