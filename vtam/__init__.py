#!/usr/bin/env python

import os
import subprocess
import sys

from vtam.utils.ArgParser import ArgParser
from vtam.CommandBlastCOI import CommandBlastCOI
from vtam.CommandTaxonomy import CommandTaxonomy
from vtam.utils.Logger import Logger
from vtam.utils.OptionManager import OptionManager
from vtam.CommandPoolMarkers import CommandPoolMarkers
from vtam.CommandTaxAssign import CommandTaxAssign
from vtam.utils.VTAMexception import VTAMexception
from vtam.utils.WopmarsRunner import WopmarsRunner

class VTAM(object):

    usage_message = """usage: vtam <command> [<args>]

        These are the VTAM commands:

   merge      Merge paired-end reads
   filter     Run VTAM filters to create an ASV table with the variant sequences as last column
   optimize   Find optimal parameters for filtering
   pool_markers   Pool overlapping markers from the ASV table into one
   taxonomy   Create the taxonomy TSV file used to create lineages 
   coi_db   Download precomputed COI Blast database 

"""

    def __init__(self, sys_argv):

        ################################################################################################################
        #
        # Parse arguments
        #
        ################################################################################################################

        self.sys_argv = sys_argv
        # AG do not use abspath for the moment. Maybe later it can be used as option
        parser = ArgParser.get_arg_parser()
        self.args = parser.parse_args(sys_argv)

        ################################################################################################################
        #
        # Add argparser attributes to optionmanager
        #
        ################################################################################################################

        option_dic = vars(self.args)
        OptionManager.instance().add_options(option_dic) # Add options to OptionManager

        # Some arguments will be passed through environmental variables
        if 'threads' in vars(self.args):
            os.environ['VTAM_THREADS'] = str(vars(self.args)['threads'])
        if 'fastqdir' in vars(self.args):
            os.environ['VTAM_FASTQ_DIR'] = vars(self.args)['fastqdir']
        if 'fastadir' in vars(self.args):
            os.environ['VTAM_FASTA_DIR'] = vars(self.args)['fastadir']

        #################################################################
        #
        # Create FilterLFNreference table
        #
        #################################################################

        if vars(self.args)['command'] in ['asv', 'optimize']:
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
                try:
                    conn.execute(filter_lfn_reference.insert(), FilterLFNreference_records)
                except:
                    pass

        ###############################################################
        #
        # Subcommands: wopfile-dependent, merge, filter, optimize
        #
        ###############################################################

        if vars(self.args)['command'] in ['merge', 'filter', 'optimize']:

            wopmars_runner = WopmarsRunner(command=vars(self.args)['command'], parameters=OptionManager.instance())
            wopmars_command = wopmars_runner.get_wopmars_command()

            ###############################################################
            #
            # Create wopmars command and implicitely wopfile
            #
            ###############################################################

            Logger.instance().info(wopmars_command)
            run_result = subprocess.run(wopmars_command, shell=True)
            sys.exit(run_result.returncode)

        ###############################################################
        #
        # Subcommand: taxassign
        #
        ###############################################################

        elif vars(self.args)['command'] == 'taxassign':
            db = OptionManager.instance()['db']
            variants_tsv = OptionManager.instance()['variants']
            variant_taxa_tsv = OptionManager.instance()['variant_taxa']
            mode = OptionManager.instance()['mode']
            taxonomy_tsv = OptionManager.instance()['taxonomy']
            blasdb_dir_path = OptionManager.instance()['blastdbdir']
            blastdbname_str = OptionManager.instance()['blastdbname']
            ltg_rule_threshold = OptionManager.instance()['ltg_rule_threshold']
            include_prop = OptionManager.instance()['include_prop']
            min_number_of_taxa = OptionManager.instance()['min_number_of_taxa']
            num_threads = OptionManager.instance()['threads']
            CommandTaxAssign.main(db=db, mode=mode, variants_tsv=variants_tsv, variant_taxa_tsv=variant_taxa_tsv, taxonomy_tsv=taxonomy_tsv,
                                  blasdb_dir_path=blasdb_dir_path, blastdbname_str=blastdbname_str,
                                  ltg_rule_threshold=ltg_rule_threshold, include_prop=include_prop,
                                  min_number_of_taxa=min_number_of_taxa, num_threads=num_threads)

        ###############################################################
        #
        # Subcommand: pool_markers
        #
        ###############################################################

        elif vars(self.args)['command'] == 'pool_markers':
            db = OptionManager.instance()['db']
            run_marker_tsv = OptionManager.instance()['runmarker']
            pooled_marker_tsv = OptionManager.instance()['pooledmarkers']
            taxonomy_tsv = OptionManager.instance()['taxonomy']
            CommandPoolMarkers.main(db=db, pooled_marker_tsv=pooled_marker_tsv, taxonomy_tsv=taxonomy_tsv,
                                    run_marker_tsv=run_marker_tsv)

        ###############################################################
        #
        # Subcommand: taxonomy
        #
        ###############################################################

        elif vars(self.args)['command'] == 'taxonomy':
            output = OptionManager.instance()['output']
            precomputed = OptionManager.instance()['precomputed']
            taxonomydb = CommandTaxonomy(output=output, precomputed=precomputed, )
            if precomputed:
                taxonomydb.download_taxonomy_tsv()
            else:
                taxonomydb.create_taxonomy_db()

        ###############################################################
        #
        # Subcommand: coi blast
        #
        ###############################################################

        elif vars(self.args)['command'] == 'coi_blast_db':
            coi_blast_db = OptionManager.instance()['coi_blast_db']
            coi_blast_db = CommandBlastCOI(coi_blast_db=coi_blast_db)
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

