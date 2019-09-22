import os
from unittest import TestCase

from wopmetabarcoding import VTAM
from wopmetabarcoding.utils.ArgParser import ArgParser


class TestParser(TestCase):

    def setUp(self):
        pass

    def test_parser_vtam_base(self):
        args_str = 'command --db test '
        #
        # Base parser
        parser = ArgParser.get_arg_parser_base()
        args = parser.parse_args(args_str.split())
        self.assertTrue(hasattr(args, 'db'))

    def test_parser_vtam_merge(self):
        args_str = 'command --db test --fastqinfo {} --fastqdir {} --fastainfo foo --fastadir foo'.format(__file__,
                os.path.dirname(__file__))
        #
        # Base parser
        parser = ArgParser.get_arg_parser_merge()
        args = parser.parse_args(args_str.split())
        self.assertTrue(hasattr(args, 'db'))

    def test_cli(self):
        args_str = 'vtam merge --db test --fastqinfo {} --fastqdir {} --fastainfo foo --fastadir foo'.format(__file__,
                                                                                                  os.path.dirname(
                                                                                                      __file__))
        VTAM(args_str.split())
