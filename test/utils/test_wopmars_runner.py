import os
from unittest import TestCase
from wopmetabarcoding.utils.ArgParser import ArgParser

from wopmetabarcoding import VTAM

from wopmetabarcoding.utils.OptionManager import OptionManager
from wopmetabarcoding.utils.WopmarsRunner import WopmarsRunner


class TestWorpmarsRunner(TestCase):

    def test_wopmars_runner_merge(self):
        #
        # Minimal merge command
        args_str = 'merge --fastqinfo {} --fastqdir {} --fastainfo foo --fastadir foo'.format(__file__,
                os.path.dirname(__file__))
        parser = ArgParser.get_arg_parser()
        args = parser.parse_args(args_str.split())
        #####################
        #
        # Add argparser attributes to optionmanager
        #
        #####################
        for k in vars(args):
            print(k,vars(args)[k])
            OptionManager.instance()[k] = vars(args)[k]
        #####################
        #
        # Test wopfile
        #
        #####################
        wopmars_runner = WopmarsRunner(tool='merge', parameters=OptionManager.instance())
        wopfile_content = wopmars_runner.get_wopfile()
        wopfile_content_bak = """rule Merge:
  tool: wopmetabarcoding.wrapper.Merge
  input:
      file:
          sample2fastq: wopmetabarcodin/test/utils/test_wopmars_runner.py
  output:
      file:
          sample2fasta: wopmetabarcodin/foo
  params:
      fastq_directory: wopmetabarcodin/test/utils
      fasta_dir: wopmetabarcodin/foo
      fastq_minovlen: 50
      fastq_maxmergelen: 300
      fastq_minmergelen: 100
      fastq_minlen: 50
      fastq_maxee: 1
      fastq_truncqual: 10
      fastq_maxns: 0
      threads: 8
      fastq_ascii: 33"""
        self.assertTrue(wopfile_content[0:95] == wopfile_content_bak[0:95])
        self.assertTrue(wopfile_content[-95:] == wopfile_content_bak[-95:])
        #####################
        #
        # Test wopmars command
        #
        #####################
        wopmars_command = wopmars_runner.get_wopmars_command()

        self.assertTrue(wopmars_command == 'wopmars -w /tmp/tmplfm8n_c7/Wopfile_merge.yml -D '
                                           'sqlite:////home/gonzalez/Software/repositories/wopmetabarcodin/db.sqlite -p '
                                           '-v --fastainfo /home/gonzalez/Software/repositories/wopmetabarcodin/foo '
                                           '--fastadir /home/gonzalez/Software/repositories/wopmetabarcodin/test/utils')


