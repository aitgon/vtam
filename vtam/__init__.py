#!/usr/bin/env python

import os
import sys

from vtam.utils.ArgParser import ArgParser
from vtam.utils.DBblastCOI import DBblastCOI
from vtam.utils.DBtaxonomy import DBtaxonomy
from vtam.utils.Logger import Logger
from vtam.utils.OptionManager import OptionManager
from vtam.utils.PoolMarkerRunner import PoolMarkerRunner
from vtam.utils.VTAMexception import VTAMexception
from vtam.utils.WopmarsRunner import WopmarsRunner

class VTAM(object):

    usage_message = """usage: vtam <command> [<args>]

        These are the VTAM commands:

   merge      Merge paired-end reads
   asv        Run VTAM from sorting reads to biosamples/replicates till making ASV table
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
        parser = ArgParser.get_arg_parser(is_abspath=False)
        self.args = parser.parse_args(sys_argv)

        #####################
        #
        # Add argparser attributes to optionmanager
        #
        #####################

        option_dic = vars(self.args)
        OptionManager.instance().add_options(option_dic) # Add options to OptionManager

        # Some arguments will be passed through environmental variables
        if 'threads' in vars(self.args):
            os.environ['VTAM_THREADS'] = str(vars(self.args)['threads'])
        if 'fastqdir' in vars(self.args):
            os.environ['VTAM_FASTQ_DIR'] = vars(self.args)['fastqdir']
        if 'fastadir' in vars(self.args):
            os.environ['VTAM_FASTA_DIR'] = vars(self.args)['fastadir']

        ###############################################################
        #
        # Subcommands: wopfile-dependent, merge, asv, optimize
        #
        ###############################################################

        if vars(self.args)['command'] in ['merge', 'asv', 'optimize']:
            wopmars_runner = WopmarsRunner(command=vars(self.args)['command'], parameters=OptionManager.instance())
            wopmars_command = wopmars_runner.get_wopmars_command()
            ###############################################################
            #
            # Create wopmars command and implicitely wopfile
            #
            ###############################################################
            Logger.instance().info(wopmars_command)
            os.system(wopmars_command)
            sys.exit(0)

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
            PoolMarkerRunner.main(db=db, pooled_marker_tsv=pooled_marker_tsv, taxonomy_tsv=taxonomy_tsv,
                                  run_marker_tsv=run_marker_tsv)

        ###############################################################
        #
        # Subcommand: taxonomy
        #
        ###############################################################

        elif vars(self.args)['command'] == 'taxonomy':
            output = OptionManager.instance()['output']
            precomputed = OptionManager.instance()['precomputed']
            taxonomydb = DBtaxonomy(output=output, precomputed=precomputed, )
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
            coi_blast_db = DBblastCOI(coi_blast_db=coi_blast_db)
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

