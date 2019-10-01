import os
from unittest import TestCase
from vtam.utils.ArgParser import ArgParser

from vtam.utils.OptionManager import OptionManager
from vtam.utils.WopmarsRunner import WopmarsRunner
from vtam.utils.PathManager import PathManager


class TestWorpmarsRunner(TestCase):

    def test_wopmars_runner_merge(self):
        #
        # Minimal merge command
        args_str = 'merge --fastqinfo {} --fastqdir {} --fastainfo foo --fastadir foo'.format(
            os.path.relpath(__file__, PathManager.get_package_path()),
            os.path.relpath(os.path.dirname(__file__), PathManager.get_package_path()))
        parser = ArgParser.get_arg_parser(abspath=False)
        args = parser.parse_args(args_str.split())
        #####################
        #
        # Add argparser attributes to optionmanager
        #
        #####################
        for k in vars(args):
            OptionManager.instance()[k] = vars(args)[k]
        ###############################################################
        #
        # Test wopfile
        #
        ###############################################################
        wopmars_runner = WopmarsRunner(command='merge', parameters=OptionManager.instance())
        wopfile_path = os.path.relpath(os.path.join(PathManager.get_package_path(), "test/output/wopfile"), PathManager.get_package_path())
        wopfile_path, wopfile_content = wopmars_runner.create_wopfile(path=wopfile_path)
        wopfile_content_bak = """rule Merge:
  tool: vtam.wrapper.Merge
  input:
      file:
          fastqinfo: test/utils/test_wopmars_runner.py
  output:
      file:
          fastainfo: foo
  params:
      fastq_dir: test/utils
      fasta_dir: foo
      fastq_minovlen: 50
      fastq_maxmergelen: 300
      fastq_minmergelen: 100
      fastq_minlen: 50
      fastq_maxee: 1
      fastq_truncqual: 10
      fastq_maxns: 0
      threads: 8
      fastq_ascii: 33
      log_verbosity: 0
      log_file: None"""
        self.assertTrue(wopfile_content == wopfile_content_bak)
        ###############################################################
        #
        # Test wopmars command
        #
        ###############################################################
        wopmars_command = wopmars_runner.get_wopmars_command()
        wopmars_runner.wopfile_path = "Wopfile_merge"
        self.assertTrue(wopmars_command == 'wopmars -w test/output/wopfile -D sqlite:///db.sqlite -p -v')


