#!/usr/bin/env python

import argparse
import sys

from wopmetabarcoding.utils.ArgParser import ArgParser
from wopmetabarcoding.utils.OptionManager import OptionManager

class VTAM(object):

    def __init__(self, sys_argv):
        self.sys_argv = sys_argv
        parser = ArgParser.get_arg_parser()
        args = parser.parse_args(sys_argv)
        # use dispatch pattern to invoke method with same name
        getattr(self, args.subcommand)()

    def merge(self):
        # now that we're inside a subcommand, ignore the first
        # TWO argvs, ie the command (vtam) and the subcommand (merge)
        parser = ArgParser.get_arg_parser()
        args = parser.parse_args(self.sys_argv[1:])
        ##################################
        #
        # Store arguments in OptionsManager
        #
        ##################################
        for k in vars(args):
            if not k is 'command':
                try:
                    OptionManager.instance()[k] = vars(args)[k][0]
                except TypeError:
                    OptionManager.instance()[k] = vars(args)[k]

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

# if __name__ == '__main__':
#     VTAM()
