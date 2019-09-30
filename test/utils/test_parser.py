from unittest import TestCase
from wopmetabarcoding.utils.ArgParser import ArgParser

import os

class TestParser(TestCase):

    def test_parser_merge(self):
        #
        # Minimal merge command
        args_str = 'merge --fastqinfo {} --fastqdir {} --fastainfo foo --fastadir foo'.format(__file__,
                os.path.dirname(__file__))
        parser = ArgParser.get_arg_parser()
        parser.parse_args(args_str.split())

    def test_parser_otu(self):
        #
        # Minimal otu command
        args_str = 'otu --fastainfo {input} --fastadir {dirname} --taxonomy {input} --outdir foo -n'.format(input=__file__,
                dirname=os.path.dirname(__file__))
        parser = ArgParser.get_arg_parser()
        parser.parse_args(args_str.split())

    def test_parser_optimize(self):
        #
        # Minimal optimize command
        args_str = 'optimize --fastainfo {} --variant_known {} --outdir foo -n'.format(__file__, __file__)
        # args_str = 'optimize --variant_known test'.format(__file__, os.path.dirname(__file__))
        parser = ArgParser.get_arg_parser()
        parser.parse_args(args_str.split())
