import os
from unittest import TestCase

from wopmetabarcoding import VTAM
from wopmetabarcoding.utils.ArgParser import ArgParser
from wopmetabarcoding.utils.OptionManager import OptionManager


class TestParser(TestCase):

    def setUp(self):
        pass

    # def test_parser_base(self):
    #     args_str = 'merge --db test '
    #     #
    #     # Base parser
    #     parser = ArgParser.get_arg_parser_base()
    #     args = parser.parse_args(args_str.split())
    #     self.assertTrue(hasattr(args, 'db'))

    def test_parser_merge(self):
        #
        # Minimal merge command
        args_str = 'merge --fastqinfo {} --fastqdir {} --fastainfo foo --fastadir foo'.format(__file__,
                os.path.dirname(__file__))
        # args_str = 'merge --db test -n'
        parser = ArgParser.get_arg_parser_merge()
        parser.parse_args(args_str.split())
        #
        # Merge command with some wopmars arguments
        args_str = 'merge --fastqinfo {} --fastqdir {} --fastainfo foo --fastadir foo -n'.format(__file__,
                os.path.dirname(__file__))
        # args_str = 'merge --db test -n'
        parser = ArgParser.get_arg_parser_merge()
        parser.parse_args(args_str.split())

    # 
    # def test_cli(self):
    #     args_str = 'vtam merge --db test --fastqinfo {} --fastqdir {} --fastainfo foo --fastadir foo'.format(__file__,
    #                                                                                               os.path.dirname(
    #                                                                                                   __file__))
    #     VTAM(args_str.split())
    # 
    # def test_1(self):
    #     OptionManager.instance()['test2'] = 'test2'
    #     print(OptionManager.instance())
    #     OptionManager.instance()['test'] = 'test3'
    #     print(OptionManager.instance())
    # 
    # def test_merge_args_to_optionsmanager(self):
    #     args_str = 'merge --db test --fastqinfo {} --fastqdir {} --fastainfo foo --fastadir foo'.format(__file__,
    #             os.path.dirname(__file__))
    #     #
    #     # Base parser
    #     parser = ArgParser.get_arg_parser_merge()
    #     args = parser.parse_args(args_str.split())
    #     for k in vars(args):
    #         if k is 'command':
    #             OptionManager.instance()[k] = vars(args)[k]
    #             continue
    #         try:
    #             OptionManager.instance()[k] = vars(args)[k][0]
    #         except TypeError:
    #             OptionManager.instance()[k] = vars(args)[k]
    #     #
    #     # assert
    #     #
    #     self.assertTrue(sorted(list(OptionManager.instance().keys())) == sorted(['command', 'dryrun', 'forceall',
    #             'log', 'params', 'targetrule', 'verbose', 'db', 'fastqinfo', 'fastainfo', 'fastqdir', 'fastadir']))
