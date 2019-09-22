#!/usr/bin/env python

import argparse
import sys

from wopmetabarcoding.utils.ArgParser import ArgParser


class VTAM(object):

    def __init__(self, sys_argv):
        self.sys_argv = sys_argv
        parser = ArgParser.get_arg_parser_base()
        args = parser.parse_args(sys_argv[1:2])
        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()

    def merge(self):
        # parser = argparse.ArgumentParser(
        #     description='Record changes to the repository')
        # # prefixing the argument with -- means it's optional
        # parser.add_argument('--amend', action='store_true')
        # now that we're inside a subcommand, ignore the first
        # TWO argvs, ie the command (git) and the subcommand (commit)
        parser = ArgParser.get_arg_parser_merge()
        # import pdb; pdb.set_trace()
        args = parser.parse_args(self.sys_argv[1:])
        ##################################
        #
        # Store in OptionsManager
        #
        ##################################
        print('Running: {}'.format(" ".join(self.sys_argv[0:])))

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
        print('Running: {}'.format(" ".join(self.sys_argv[0:])))

if __name__ == '__main__':
    VTAM()
