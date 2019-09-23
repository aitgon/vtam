import os
from unittest import TestCase
from wopmetabarcoding.utils.ArgParser import ArgParser

from wopmetabarcoding import VTAM

from wopmetabarcoding.utils.OptionManager import OptionManager


class TestParser(TestCase):

    def test_parser_merge(self):
        #
        # Minimal merge command
        args_str = 'merge --fastqinfo {} --fastqdir {} --fastainfo foo --fastadir foo'.format(__file__,
                os.path.dirname(__file__))
        parser = ArgParser.get_arg_parser()
        parser.parse_args(args_str.split())

        # Merge command with some wopmars arguments
        args_str = 'merge --fastqinfo {} --fastqdir {} --fastainfo foo --fastadir foo -n'.format(__file__,
                os.path.dirname(__file__))
        parser = ArgParser.get_arg_parser()
        parser.parse_args(args_str.split())

        # Run VTAM and store OptionManager
        VTAM(args_str.split())
        self.assertTrue(sorted(list(OptionManager.instance().keys())) == ['db', 'dryrun', 'fastadir', 'fastainfo', 'fastqdir',
                                                          'fastqinfo', 'forceall', 'log', 'params', 'targetrule',
                                                          'verbose'])


    def test_parser_otu(self):
        #
        # Minimal otu command
        args_str = 'otu'.format(__file__,
                os.path.dirname(__file__))
        parser = ArgParser.get_arg_parser()
        parser.parse_args(args_str.split())
        #
        # Typical otu command
        args_str = 'otu --outdir test2'.format(__file__,
                os.path.dirname(__file__))
        parser = ArgParser.get_arg_parser()
        parser.parse_args(args_str.split())

    def test_parser_optimize(self):
        #
        # Minimal optimize command
        args_str = 'optimize --variant_known test'.format(__file__,
                os.path.dirname(__file__))
        parser = ArgParser.get_arg_parser()
        parser.parse_args(args_str.split())
        #
        # Typical optimize command
        args_str = 'optimize --variant_known test --outdir test2'.format(__file__,
                os.path.dirname(__file__))
        parser = ArgParser.get_arg_parser()
        parser.parse_args(args_str.split())
