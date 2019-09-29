#!/usr/bin/env python

import sys

import os

from wopmetabarcoding.utils.ArgParser import ArgParser
from wopmetabarcoding.utils.Logger import Logger
from wopmetabarcoding.utils.OptionManager import OptionManager
from wopmetabarcoding.utils.VTAMexception import VTAMexception
from wopmetabarcoding.utils.WopmarsRunner import WopmarsRunner


class VTAM(object):

    def __init__(self, sys_argv):
        self.sys_argv = sys_argv
        parser = ArgParser.get_arg_parser(abspath=True)
        self.args = parser.parse_args(sys_argv)
        #####################
        #
        # Add argparser attributes to optionmanager
        #
        #####################
        for k in vars(self.args):
            OptionManager.instance()[k] = vars(self.args)[k]
        try:
            wopmars_runner = WopmarsRunner(command=vars(self.args)['command'], parameters=OptionManager.instance())
            wopmars_command = wopmars_runner.get_wopmars_command()
            ###############################################################
            #
            # Create wopmars command and implicitely wopfile
            #
            ###############################################################
            Logger.instance().info(wopmars_command)
            import pdb; pdb.set_trace()
            os.system(wopmars_command)
            sys.exit(0)
        except KeyError:
            print(VTAMexception(message="""usage: vtam <command> [<args>]

        These are the VTAM commands:

   merge      Merge FASTQ files
   otu        Carry out the whole pipeline, including sort and count reads, filter variants, tax assign and create OTU table
   optimize   Show different variant characteristics to help select filter parameters
"""))
