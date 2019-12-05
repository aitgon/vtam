#!/usr/bin/env python

import os
import sys

from vtam.utils.ArgParser import ArgParser
from vtam.utils.Logger import Logger
from vtam.utils.OptionManager import OptionManager
from vtam.utils.VTAMexception import VTAMexception
from vtam.utils.WopmarsRunner import WopmarsRunner


class VTAM(object):

    usage_message = """usage: vtam <command> [<args>]

        These are the VTAM commands:

   merge      Merge FASTQ files
   asv        Carry out the whole pipeline, including sort and count reads, filter variants, tax assign and create ASV table
   optimize   Show different variant characteristics to help select filter parameters
"""

    def __init__(self, sys_argv):

        ################################################################################################################
        #
        # Parse arguments
        #
        ################################################################################################################

        self.sys_argv = sys_argv
        parser = ArgParser.get_arg_parser(is_abspath=True)
        self.args = parser.parse_args(sys_argv)

        #####################
        #
        # Add argparser attributes to optionmanager
        #
        #####################

        for k in vars(self.args):
            OptionManager.instance()[k] = vars(self.args)[k]
        # Some arguments will be passed through environmental variables
        if 'threads' in vars(self.args):
            os.environ['VTAM_THREADS'] = str(vars(self.args)['threads'])
        if 'fastqdir' in vars(self.args):
            os.environ['VTAM_FASTQ_DIR'] = vars(self.args)['fastqdir']
        if 'fastadir' in vars(self.args):
            os.environ['VTAM_FASTA_DIR'] = vars(self.args)['fastadir']
        try:
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
        except KeyError:
            Logger.instance().error(VTAMexception(message=VTAM.usage_message))


def main():
    if sys.argv[1:] == []:  # No arguments
        Logger.instance().error(VTAMexception(message=VTAM.usage_message))
        sys.exit(1)
    VTAM(sys.argv[1:])

