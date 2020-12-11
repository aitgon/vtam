__version__ = "0.1.21"

import os
import pathlib
import sys
import yaml

from vtam import CommandExample
from vtam.CommandBlastCOI import CommandBlastCOI
from vtam.CommandExample import CommandExample
from vtam.CommandFilterOptimize import CommandFilterOptimize
from vtam.CommandMerge import CommandMerge
from vtam.CommandPoolRunMarkers import CommandPoolRunMarkers
from vtam.CommandSortReads import CommandSortReads
from vtam.CommandTaxAssign import CommandTaxAssign
from vtam.CommandTaxonomy import CommandTaxonomy
from vtam.utils.ArgParser import ArgParser
from vtam.utils.FileCutoffSpecific import FileCutoffSpecific
from vtam.utils.Logger import Logger
from vtam.utils.Logger import LoggerArguments
from vtam.utils.FileParams import FileParams
from vtam.utils.PathManager import PathManager
from vtam.utils.VTAMexception import VTAMexception
from vtam.utils.RunnerWopmars import RunnerWopmars
from vtam.utils.constants import FilterLFNreference_records


class VTAM(object):

    def __init__(self, sys_argv):

        ############################################################################################
        #
        # Parse arguments
        #
        ############################################################################################

        self.sys_argv = sys_argv
        # AG do not use abspath for the moment. Maybe later it can be used as
        # option
        parser = ArgParser.get_main_arg_parser()
        self.args = parser.parse_args(sys_argv)

        arg_parser_dic = vars(self.args)

        ############################################################################################
        #
        # If non-specified, initiate params.yml
        #
        ############################################################################################

        if 'params' in arg_parser_dic and arg_parser_dic['params'] is None:
            params_yml = os.path.join(PathManager.instance().get_configdir(), "params.yml")
            if not os.path.isfile(params_yml):
                pathlib.Path(params_yml).touch(exist_ok=False)
            arg_parser_dic['params'] = params_yml

        ############################################################################################
        #
        # Parse log arguments
        #
        ############################################################################################

        if 'log_verbosity' in arg_parser_dic:
            (LoggerArguments.instance()).update({'log_verbosity': arg_parser_dic['log_verbosity']})
            os.environ['VTAM_LOG_VERBOSITY'] = str(
                arg_parser_dic['log_verbosity'])

        if 'log' in arg_parser_dic:
            (LoggerArguments.instance()).update({'log': arg_parser_dic['log']})
            os.environ['VTAM_LOG_FILE'] = str(arg_parser_dic['log'])

        #######################################################################
        #
        # Set arguments, logger
        #
        #######################################################################

        # Some arguments will be passed through environmental variables
        if 'threads' in arg_parser_dic:
            os.environ['VTAM_THREADS'] = str(arg_parser_dic['threads'])

        ############################################################################################
        #
        # Subcommands: wopfile-dependent, filter, optimize
        #
        ############################################################################################

        if arg_parser_dic['command'] in ['filter', 'optimize']:

            if arg_parser_dic['command'] in ['filter']:

                ####################################################################################
                #
                # Verify coherence of --lfn_variant_replicate and params arguments
                #
                ####################################################################################

                with open(arg_parser_dic['params']) as fin:
                    # The FullLoader parameter handles the conversion from YAML
                    # scalar values to Python the dictionary format
                    params_dic = yaml.load(fin, Loader=yaml.SafeLoader) or {}

                    if arg_parser_dic['lfn_variant_replicate']:
                        if 'lfn_variant_cutoff' in params_dic:
                            Logger.instance().error(VTAMexception(
                                'The parameter "lfn_variant_cutoff" in the parameter file "{}" is incompatible with'
                                ' the --lfn_variant_replicate argument.'.format(arg_parser_dic['params'])))
                            sys.exit(1)

                    else:
                        if 'lfn_variant_replicate_cutoff' in params_dic:
                            Logger.instance().error(VTAMexception(
                                'The parameter "lfn_variant_replicate_cutoff" in the parameter file "{}" needs'
                                ' the --lfn_variant_replicate argument.'.format(arg_parser_dic['params'])))
                            sys.exit(1)

                ####################################################################################
                #
                # Verify coherence of --lfn_variant_replicate and cutoff_specific argument
                #
                ####################################################################################

                if not (arg_parser_dic['cutoff_specific'] is None):  # cutoff specific argument

                    if arg_parser_dic['lfn_variant_replicate']:  # lfn_variant_replicate

                        # cutoff_specific for lfn_variant
                        if not FileCutoffSpecific(arg_parser_dic['cutoff_specific']).is_compatible_lfn_variant_replicate():
                            Logger.instance().error('The --lfn_variant_replicate argument is incompatible with the cutoff_specific file {}.'.format(
                                    arg_parser_dic['cutoff_specific']))
                            sys.exit(1)

                    else: # lfn_variant

                        # cutoff_specific for lfn_variant_replicate
                        if FileCutoffSpecific(arg_parser_dic['cutoff_specific']).is_compatible_lfn_variant_replicate():

                            Logger.instance().error('The cutoff_specific file {} requires the --lfn_variant_replicate argument.'.format(
                                    arg_parser_dic['cutoff_specific']))
                            sys.exit(1)

                ############################################################################################
                #
                # If non-specified, initiate cutoff specific
                #
                ############################################################################################

                if arg_parser_dic['cutoff_specific'] is None:
                    cutoff_specific_tsv = os.path.join(PathManager.instance().get_configdir(),
                                                       "cutoff_specific.tsv")
                    if not os.path.isfile(cutoff_specific_tsv):
                        pathlib.Path(cutoff_specific_tsv).touch(exist_ok=False)
                    arg_parser_dic['cutoff_specific'] = cutoff_specific_tsv

            CommandFilterOptimize.main(arg_parser_dic=arg_parser_dic)

        ############################################################################################
        #
        # Subcommand: example
        #
        ############################################################################################

        elif arg_parser_dic['command'] == 'example':
            outdir = arg_parser_dic['outdir']
            CommandExample.main(outdir=outdir)

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
            sorteddir = arg_parser_dic['sorteddir']
            num_threads = arg_parser_dic['threads']
            params = arg_parser_dic['params']
            CommandSortReads.main(fastainfo=fastainfo, fastadir=fastadir, params=params,
                                  num_threads=num_threads, sorteddir=sorteddir)

        ############################################################################################
        #
        # Subcommand: taxassign
        #
        ############################################################################################

        elif arg_parser_dic['command'] == 'taxassign':
            db = arg_parser_dic['db']
            asvtable_tsv = arg_parser_dic['asvtable']
            output = arg_parser_dic['output']
            mode = arg_parser_dic['mode']
            taxonomy_tsv = arg_parser_dic['taxonomy']
            blasdb_dir_path = arg_parser_dic['blastdbdir']
            blastdbname_str = arg_parser_dic['blastdbname']
            num_threads = arg_parser_dic['threads']
            params = arg_parser_dic['params']
            CommandTaxAssign.main(db=db, mode=mode, asvtable_tsv=asvtable_tsv, output=output,
                                  taxonomy_tsv=taxonomy_tsv, blastdb_dir_path=blasdb_dir_path,
                                  blastdbname_str=blastdbname_str, params=params, num_threads=num_threads)

        ############################################################################################
        #
        # Subcommand: pool
        #
        ############################################################################################

        elif arg_parser_dic['command'] == 'pool':
            db = arg_parser_dic['db']
            readcounts = arg_parser_dic['readcounts']
            run_marker_tsv = arg_parser_dic['runmarker']
            pooled_marker_tsv = arg_parser_dic['asvtable']
            params = arg_parser_dic['params']
            CommandPoolRunMarkers.main(db=db, pooled_marker_tsv=pooled_marker_tsv,
                run_marker_tsv=run_marker_tsv, params=params, readcounts=readcounts)

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
            blastdbdir = arg_parser_dic['blastdbdir']
            blastdbname = arg_parser_dic['blastdbname']
            coi_blast_db = CommandBlastCOI(blastdbname=blastdbname)
            coi_blast_db.download(blastdbdir=blastdbdir)

        ############################################################################################
        #
        # Else: run_name usage message
        #
        ############################################################################################

        else:
            self.args = parser.parse_args(['--help'])  # if command unknown print help



def main():

    if not sys.argv[1:]:  # if not arguments, print help
        VTAM(['--help'])
    VTAM(sys.argv[1:])
