import os
import pathlib
import shlex
import unittest
import yaml

from vtam.utils.ArgParser import ArgParser
from vtam.utils.PathManager import PathManager
from vtam.utils.WopmarsRunner import WopmarsRunner


class TestWorpmarsRunnerFilter(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        package_path = PathManager.get_package_path()

        foopaths = {}
        foopaths['foofile'] = os.path.relpath(__file__, PathManager.get_package_path())
        foopaths['foodir'] = os.path.relpath(os.path.dirname(__file__), PathManager.get_package_path())
        foopaths['outdir'] = os.path.relpath(os.path.join(PathManager.get_test_path(), 'output'),
            PathManager.get_package_path())
        foopaths['blastdb'] = os.path.relpath(os.path.join(PathManager.get_test_path(), 'test_files', 'blastdb'),
            PathManager.get_package_path())
        foopaths['readinfo_tsv'] = os.path.relpath(os.path.join(package_path, "doc/data/readinfo_mfzr.tsv"),
            PathManager.get_package_path())
        cls.foopaths = foopaths

        cls.minseqlength_value_32 = 32
        cls.minseqlength_value_40 = 40
        cls.lfn_variant_replicate_cutoff = 0.002

    def setUp(self):

        self.tempdir = PathManager.instance().get_tempdir()
        pathlib.Path(self.tempdir).mkdir(parents=True, exist_ok=True)

    def test_wopmars_runner_filter(self):

        args_str = 'filter --readinfo {readinfo_tsv} --readdir {foodir} --asvtable asvtableoutput.tsv'.format(**self.foopaths)
        args = ArgParser.get_main_arg_parser().parse_args(shlex.split(args_str))

        wopmars_runner = WopmarsRunner(command='filter', cli_args_dic=vars(args))
        wopfile_path, wopfile_content = wopmars_runner.create_wopfile()

        with open(os.path.join(os.path.dirname(__file__), "wopfile_filter.yml")) as fin:
            wopfile_content_bak = fin.read()

        self.assertTrue(wopfile_content == wopfile_content_bak.strip())
        self.assertTrue('lfn_variant_cutoff' in yaml.load(wopfile_content, Loader=yaml.SafeLoader)['rule FilterLFN']['params'])

    def test_wopmars_runner_filter_lfn_variant_replicate(self):

        args_str = 'filter --readinfo {readinfo_tsv} --readdir {foodir} --asvtable asvtableoutput.tsv --lfn_variant_replicate'.format(
            **self.foopaths)
        args = ArgParser.get_main_arg_parser().parse_args(shlex.split(args_str))

        wopmars_runner = WopmarsRunner(command='filter', cli_args_dic=vars(args))
        wopfile_path, wopfile_content = wopmars_runner.create_wopfile()

        self.assertFalse('lfn_variant_cutoff' in yaml.load(wopfile_content, Loader=yaml.SafeLoader)['rule FilterLFN']['params'])
        self.assertTrue('lfn_variant_replicate_cutoff' in yaml.load(wopfile_content, Loader=yaml.SafeLoader)['rule FilterLFN']['params'])

    def test_wopmars_runner_filter_with_cutoff_specific(self):

        args_str = 'filter --readinfo {readinfo_tsv} --readdir {foodir} --asvtable asvtableoutput.tsv' \
                   ' --cutoff_specific {foofile}'.format(**self.foopaths)
        args = ArgParser.get_main_arg_parser().parse_args(args_str.split())

        wopmars_runner = WopmarsRunner(command='filter', cli_args_dic=vars(args))
        wopfile_path = os.path.relpath(os.path.join(PathManager.get_package_path(), "tests/output/wopfile"), PathManager.get_package_path())
        wopfile_path, wopfile_content = wopmars_runner.create_wopfile(path=wopfile_path)

        self.assertTrue(yaml.load(wopfile_content, Loader=yaml.SafeLoader)['rule FilterLFN']['input']['file']['cutoff_specific'] == self.foopaths['foofile'])
