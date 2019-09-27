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
        # import pdb; pdb.set_trace()

        # Merge command with some wopmars arguments
        args_str = 'merge --fastqinfo {} --fastqdir {} --fastainfo foo --fastadir foo -n'.format(__file__,
                os.path.dirname(__file__))
        parser = ArgParser.get_arg_parser()
        parser.parse_args(args_str.split())
        #
        # Run VTAM and store OptionManager
        # import pdb; pdb.set_trace()
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

    def test_wopfile_merge(self):
        """rule Merge:
  tool: wopmetabarcoding.wrapper.Merge
  input:
      file:
          sample2fastq: /home/gonzalez/Software/repositories/wopmetabarcodin/test/test_parser.py
  output:
      file:
          sample2fasta: /home/gonzalez/Software/repositories/wopmetabarcodin/foo
  params:
      fastq_directory: /home/gonzalez/Software/repositories/wopmetabarcodin/test
      fasta_dir: /home/gonzalez/Software/repositories/wopmetabarcodin/foo
      fastq_minovlen: 50
      fastq_maxmergelen: 300
      fastq_minmergelen: 100
      fastq_minlen: 50
      fastq_maxee: 1
      fastq_truncqual: 10
      fastq_maxns: 0
      threads: 8
      fastq_ascii: 33
      """
