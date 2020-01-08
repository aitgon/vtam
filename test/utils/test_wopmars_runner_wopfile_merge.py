import os
import pathlib
from unittest import TestCase
from vtam.utils.ArgParser import ArgParser

from vtam.utils.OptionManager import OptionManager
from vtam.utils.WopmarsRunner import WopmarsRunner
from vtam.utils.PathManager import PathManager


class TestWorpmarsRunnerMerge(TestCase):

    def test_wopmars_runner_merge(self):
        #
        # Minimal merge command
        args_str = 'merge --fastqinfo {} --fastqdir {} --fastainfo foo --fastadir foo'.format(
            os.path.relpath(__file__, PathManager.get_package_path()),
            os.path.relpath(os.path.dirname(__file__), PathManager.get_package_path()))
        parser = ArgParser.get_arg_parser(is_abspath=False)
        args = parser.parse_args(args_str.split())

        #####################
        #
        # Add argparser attributes to optionmanager
        #
        #####################

        option_dic = vars(args)
        OptionManager.instance().add_options(option_dic) # Add options to OptionManager

        ###############################################################
        #
        # Test wopfile
        #
        ###############################################################
        wopmars_runner = WopmarsRunner(command='merge', parameters=OptionManager.instance())
        wopfile_path = os.path.relpath(os.path.join(PathManager.get_package_path(), "test/output/wopfile"),
                                    PathManager.get_package_path())
        wopfile_path, wopfile_content = wopmars_runner.create_wopfile(path=wopfile_path)
        wopfile_content_bak = """rule Merge:
    tool: vtam.wrapper.Merge
    input:
        file:
            fastqinfo: test/utils/test_wopmars_runner_wopfile_merge.py
    output:
        file:
            fastainfo: foo
    params:
        fastq_minovlen: 50
        fastq_maxmergelen: 500
        fastq_minmergelen: 100
        fastq_minlen: 50
        fastq_maxee: 1
        fastq_truncqual: 10
        fastq_maxns: 0
        fastq_ascii: 33"""
        self.assertTrue(wopfile_content == wopfile_content_bak)
        ###############################################################
        #
        # Test wopmars command
        #
        ###############################################################
        wopmars_command = wopmars_runner.get_wopmars_command()
        wopmars_runner.wopfile_path = "Wopfile_merge"
        self.assertTrue(wopmars_command == 'wopmars -w test/output/wopfile -D sqlite:///db.sqlite ')


    def test_wopmars_runner_merge_change_params(self):

        #####################
        #
        # Params yml
        #
        #####################

        tempdir = PathManager.instance().get_tempdir()
        pathlib.Path(tempdir).mkdir(parents=True, exist_ok=True)
        params_yml_str = "fastq_minovlen: 100"
        params_yml_path = os.path.join(tempdir, "params.yml")
        with open(params_yml_path, "w") as fout:
            fout.write(params_yml_str)
        #
        args_str = 'merge --fastqinfo {} --fastqdir {} --fastainfo foo --fastadir foo -vv --params {}'.format(
            os.path.relpath(__file__, PathManager.get_package_path()),
            os.path.relpath(os.path.dirname(__file__), PathManager.get_package_path()),
            params_yml_path)
        parser = ArgParser.get_arg_parser(is_abspath=False)
        args = parser.parse_args(args_str.split())
        option_dic = vars(args) # Dictionnary with options
        #####################
        #
        # Add argparser attributes to optionmanager
        #
        #####################
        OptionManager.instance().add_options(option_dic) # Add options to OptionManager
        ###############################################################
        #
        # Test wopfile
        #
        ###############################################################
        wopmars_runner = WopmarsRunner(command='merge', parameters=OptionManager.instance())
        wopfile_path = os.path.relpath(os.path.join(PathManager.get_package_path(), "test/output/wopfile"),
                                    PathManager.get_package_path())
        wopfile_path, wopfile_content = wopmars_runner.create_wopfile(path=wopfile_path)
        wopfile_content_bak = """rule Merge:
    tool: vtam.wrapper.Merge
    input:
        file:
            fastqinfo: test/utils/test_wopmars_runner_wopfile_merge.py
    output:
        file:
            fastainfo: foo
    params:
        fastq_minovlen: 100
        fastq_maxmergelen: 500
        fastq_minmergelen: 100
        fastq_minlen: 50
        fastq_maxee: 1
        fastq_truncqual: 10
        fastq_maxns: 0
        fastq_ascii: 33"""
        self.assertTrue(wopfile_content == wopfile_content_bak)

