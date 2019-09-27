#!/usr/bin/env python

import argparse
import sys

import os

from wopmetabarcoding.utils.ArgParser import ArgParser
from wopmetabarcoding.utils.OptionManager import OptionManager
from wopmetabarcoding.utils.VTAMexception import VTAMexception
from wopmetabarcoding.utils.WopmarsRunner import WopmarsRunner


class VTAM(object):

    def __init__(self, sys_argv):
        self.sys_argv = sys_argv
        parser = ArgParser.get_arg_parser()
        args = parser.parse_args(sys_argv)
        # use dispatch pattern to invoke method with same name
        try:
            getattr(self, vars(args)['command'])()
        except AttributeError:
            print(VTAMexception(message="""usage: vtam <command> [<args>]
        
        These are the VTAM commands:

   merge      Merge FASTQ files
   otu        Carry out the whole pipeline, including sort and count reads, filter variants, tax assign and create OTU table
   optimize   Show different variant characteristics to help select filter parameters
"""))

    def merge(self):
        # now that we're inside a subcommand, ignore the first
        # argvs, ie the command (vtam) and the subcommand (merge)
        parser = ArgParser.get_arg_parser()
        args = parser.parse_args(self.sys_argv)
        #####################
        #
        # Add argparser attributes to optionmanager
        #
        #####################
        for k in vars(args):
            OptionManager.instance()[k] = vars(args)[k]
        ###############################################################
        #
        # Create wopmars command and implicitely wopfile
        #
        ###############################################################
        wopmars_runner = WopmarsRunner(subcommand='merge', parameters=OptionManager.instance())
        wopmars_command = wopmars_runner.get_wopmars_command()
        ###############################################################
        #
        # Create wopmars command and implicitely wopfile
        #
        ###############################################################
        os.system(wopmars_command)
        sys.exit(0)

    def otu(self):
        parser = argparse.ArgumentParser(
            description='Download objects and refs from another repository')
        # NOT prefixing the argument with -- means it's not optional
        parser.add_argument('repository')
        args = parser.parse_args(sys.argv[2:])
        ##################################
        #
        # Store in OptionsManager
        #
        ##################################
        for k in vars(args):
            if not k is 'command':
                try:
                    OptionManager.instance()[k] = vars(args)[k][0]
                except TypeError:
                    OptionManager.instance()[k] = vars(args)[k]
