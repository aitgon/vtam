#!/usr/bin/env python

import os
import subprocess
import sys

import sqlalchemy

from vtam.CommandMerge import CommandMerge
from vtam.CommandSortReads import CommandSortReads
from vtam.utils.ArgParser import ArgParser
from vtam.CommandBlastCOI import CommandBlastCOI
from vtam.CommandTaxonomy import CommandTaxonomy
from vtam.utils.Logger import Logger
from vtam.CommandPoolRunMarkers import CommandPoolRunMarkers
from vtam.CommandTaxAssign import CommandTaxAssign
from vtam.utils.VTAMexception import VTAMexception
from vtam.utils.WopmarsRunner import WopmarsRunner
from vtam.utils.Logger import LoggerArguments

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

        ################################################################################################################
        #
        # Parse arguments
        #
        ################################################################################################################

        self.sys_argv = sys_argv
        # AG do not use abspath for the moment. Maybe later it can be used as option
        parser = ArgParser.get_main_arg_parser()
        self.args = parser.parse_args(sys_argv)

        arg_parser_dic = vars(self.args)

        ################################################################################################################
        #
        # Parse log arguments
        #
        ################################################################################################################

        (LoggerArguments.instance()).update({'log_verbosity': arg_parser_dic['log_verbosity'],
                                             'log_file': arg_parser_dic['log_file']})
        if 'log_verbosity' in arg_parser_dic:
            os.environ['VTAM_LOG_VERBOSITY'] = str(arg_parser_dic['log_verbosity'])
        if 'log_file' in arg_parser_dic:
            os.environ['VTAM_LOG_FILE'] = str(arg_parser_dic['log_file'])

        ################################################################################################################
        #
        # Set arguments, logger
        #
        ################################################################################################################

        # Some arguments will be passed through environmental variables
        if 'threads' in vars(self.args):
            os.environ['VTAM_THREADS'] = str(vars(self.args)['threads'])

        ###############################################################
        #
        # Subcommands: wopfile-dependent, filter, optimize
        #
        ###############################################################

        if vars(self.args)['command'] in ['filter', 'optimize']:

            ############################################################################################################
            #
            # Filter and optimize CLI arguments and numerical parameters
            #
            ############################################################################################################

            ############################################################################################################
            #
            # Create FilterLFNreference table and fill it
            #
            ############################################################################################################

            from sqlalchemy import create_engine
            from sqlalchemy import Table, Column, Integer, String, MetaData
            from vtam.utils.constants import FilterLFNreference_records
            engine = create_engine('sqlite:///{}'.format(str(vars(self.args)['db'])), echo=False)
            meta = MetaData()
            filter_lfn_reference = Table(
                'FilterLFNreference', meta,
                Column('filter_id', Integer, primary_key=True),
                Column('filter_name', String),
            )
            meta.create_all(engine)

            with engine.connect() as conn:
                for filter_rec in FilterLFNreference_records:
                    filter_name = filter_rec['filter_name']
                    select_row = conn.execute(sqlalchemy.select([filter_lfn_reference.c.filter_id])
                                              .where(filter_lfn_reference.c.filter_name == filter_name)).first()
                    if select_row is None:  # variant_sequence IS NOT in the database, so INSERT it
                        conn.execute(filter_lfn_reference.insert().values(**filter_rec))

            wopmars_runner = WopmarsRunner(command=vars(self.args)['command'], cli_args_dic=arg_parser_dic)
            wopmars_command = wopmars_runner.get_wopmars_command()

            ############################################################################################################
            #
            # Create wopmars command and implicitely wopfile
            #
            ############################################################################################################

            # Some arguments will be passed through environmental variables
            if 'threads' in vars(self.args):
                os.environ['VTAM_THREADS'] = str(vars(self.args)['threads'])
            Logger.instance().info(wopmars_command)
            run_result = subprocess.run(wopmars_command, shell=True)
            sys.exit(run_result.returncode)

        ###############################################################
        #
        # Subcommand: merge
        #
        ###############################################################

        elif vars(self.args)['command'] == 'merge':
            fastqinfo = arg_parser_dic['fastqinfo_tsv_path']
            fastqdir = arg_parser_dic['fastqdir']
            fastainfo = arg_parser_dic['fastainfo']
            fastadir = arg_parser_dic['fastadir']
            num_threads = arg_parser_dic['threads']
            params = arg_parser_dic['params']
            CommandMerge.main(fastqinfo=fastqinfo, fastqdir=fastqdir, fastainfo=fastainfo, fastadir=fastadir,
                              params=params, num_threads=num_threads)

        ###############################################################
        #
        # Subcommand: sortreads
        #
        ###############################################################

        elif vars(self.args)['command'] == 'sortreads':
            fastadir = arg_parser_dic['fastadir']
            fastainfo = arg_parser_dic['fastainfo']
            outdir = arg_parser_dic['outdir']
            num_threads = arg_parser_dic['threads']
            params = arg_parser_dic['params']
            CommandSortReads.main(fastainfo=fastainfo, fastadir=fastadir, params=params, num_threads=num_threads, outdir=outdir)

        ###############################################################
        #
        # Subcommand: taxassign
        #
        ###############################################################

        elif vars(self.args)['command'] == 'taxassign':
            db = arg_parser_dic['db']
            variants_tsv = arg_parser_dic['variants']
            output = arg_parser_dic['output']
            mode = arg_parser_dic['mode']
            taxonomy_tsv = arg_parser_dic['taxonomy']
            blasdb_dir_path = arg_parser_dic['blastdbdir']
            blastdbname_str = arg_parser_dic['blastdbname']
            num_threads = arg_parser_dic['threads']
            params = arg_parser_dic['params']
            CommandTaxAssign.main(db=db, mode=mode, variants_tsv=variants_tsv, output=output, taxonomy_tsv=taxonomy_tsv,
                                  blasdb_dir_path=blasdb_dir_path, blastdbname_str=blastdbname_str, params=params,
                                  num_threads=num_threads)

        ###############################################################
        #
        # Subcommand: pool
        #
        ###############################################################

        elif vars(self.args)['command'] == 'pool':
            db = arg_parser_dic['db']
            run_marker_tsv = arg_parser_dic['runmarker']
            pooled_marker_tsv = arg_parser_dic['output']
            CommandPoolRunMarkers.main(db=db, pooled_marker_tsv=pooled_marker_tsv, run_marker_tsv=run_marker_tsv)

        ###############################################################
        #
        # Subcommand: taxonomy
        #
        ###############################################################

        elif vars(self.args)['command'] == 'taxonomy':
            taxonomy_tsv = arg_parser_dic['output']
            precomputed = arg_parser_dic['precomputed']
            taxonomy = CommandTaxonomy(taxonomy_tsv=taxonomy_tsv)
            taxonomy.main()

        ###############################################################
        #
        # Subcommand: coi blast
        #
        ###############################################################

        elif vars(self.args)['command'] == 'coi_blast_db':
            coi_blast_db_dir = arg_parser_dic['coi_blast_db_dir']
            coi_blast_db = CommandBlastCOI(coi_blast_db_dir=coi_blast_db_dir)
            coi_blast_db.download()

        ###############################################################
        #
        # Else: run usage message
        #
        ###############################################################

        else:
            Logger.instance().error(VTAMexception(message=VTAM.usage_message))


def main():
    if sys.argv[1:] == []:  # No arguments
        Logger.instance().error(VTAMexception(message=VTAM.usage_message))
        sys.exit(1)
    VTAM(sys.argv[1:])

