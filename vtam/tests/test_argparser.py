import os
import pathlib
import shutil
import unittest

from vtam.utils.ArgParser import ArgParser
from vtam.utils.PathManager import PathManager


class TestArgParser(unittest.TestCase):

    def setUp(self):
        self.parser = ArgParser.get_main_arg_parser()

        package_path = PathManager.get_package_path()
        test_path = PathManager.get_test_path()
        self.test_path = test_path
        outdir_path = os.path.join(test_path, "outdir")

        self.foopaths = {}
        self.foopaths['filedoesnotexist'] = "filedoesnotexist"
        self.foopaths['dirdoesnotexist'] = "dirdoesnotexist"
        self.foopaths['fileisempty'] = os.path.join(test_path, "test_files", "emptyfile")
        self.foopaths['filenottsv'] = __file__
        self.foopaths['sortedinfo_tsv'] = os.path.join("../data/example/sortedinfo_mfzr.tsv")
        self.foopaths['params_yml'] = os.path.join(package_path, "data/example", "params_mfzr.yml")
        self.foopaths['params_wrong_yml'] = os.path.join(test_path, "test_params_file", "params_wrong.yml")
        self.foopaths['known_occurrences'] = os.path.join(package_path, "data/example", "known_occurrences.tsv")
        self.foopaths['asvtable_tsv'] = os.path.join(
            test_path, "test_files_dryad.f40v5_small", "run1_mfzr_zfzr", "asvtable_default.tsv")
        self.foopaths['runmarker_tsv'] = os.path.join(package_path, "data/example", "pool_run_marker.tsv")

        self.foopaths['taxonomy_tsv'] = os.path.join(PathManager.get_test_path(),
            "test_files_dryad.f40v5_small", "taxonomy.tsv")

        self.foopaths['foodir'] = package_path
        self.foopaths['sorteddir'] = outdir_path
        self.foopaths['emptydir'] = os.path.join(outdir_path, 'emptydir')
        pathlib.Path(
            os.path.join(
                self.foopaths['emptydir'])).mkdir(
            parents=True,
            exist_ok=True)
        self.foopaths['blastdb'] = os.path.relpath(
            os.path.join(
                PathManager.get_test_path(),
                'test_files',
                'blastdb'),
            PathManager.get_package_path())

    def test_arg_parser_params(self):

        args = "filter --sortedinfo {sortedinfo_tsv} --sorteddir {foodir} --asvtable asvtable.tsv --params {params_yml}".format(
            **self.foopaths).split()
        self.assertTrue(self.parser.parse_args(args), 0)

        args = "filter --sortedinfo {sortedinfo_tsv} --sorteddir {foodir} --asvtable asvtable.tsv --params {params_wrong_yml}".format(
            **self.foopaths).split()
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args)

    def test_arg_parser_filter(self):

        # Ok
        args = "filter --sortedinfo {sortedinfo_tsv} --sorteddir {foodir} --asvtable asvtable.tsv".format(
            **self.foopaths).split()
        self.assertTrue(self.parser.parse_args(args), 0)

        #######################################################################
        #
        # raises SystemExit
        #
        #######################################################################

        args = ["filter"]
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args)

        args = "filter --sortedinfo {filedoesnotexist} --sorteddir {foodir} --sorteddir {foodir}".format(
            **self.foopaths).split()
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args)

        args = "filter --sortedinfo {filenottsv} --sorteddir {foodir} --sorteddir {foodir}".format(
            **self.foopaths).split()
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args)

        args = "filter --sortedinfo {sortedinfo_tsv} --sorteddir {emptydir} --sorteddir {foodir}".format(
            **self.foopaths).split()
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args)

        args = "filter --sortedinfo {sortedinfo_tsv} --sorteddir {dirdoesnotexist} --sorteddir {foodir}".format(
            **self.foopaths).split()
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args)

    def test_arg_parser_optimize(self):

        # Ok
        args = "optimize --known_occurrences {known_occurrences} --sortedinfo {sortedinfo_tsv} --sorteddir {foodir} " \
               "--outdir {foodir}".format(**self.foopaths).split()
        self.assertTrue(self.parser.parse_args(args), 0)

        #######################################################################
        #
        # raises SystemExit
        #
        #######################################################################

        args = ["optimize"]
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args)

        args = "optimize --known_occurrences {filedoesnotexist} --sortedinfo {sortedinfo_tsv} --sorteddir {foodir} " \
               "--sorteddir {foodir}".format(**self.foopaths).split()
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args)

        args = "optimize --known_occurrences {filenottsv} --sortedinfo {sortedinfo_tsv} --sorteddir {foodir} " \
               "--sorteddir {foodir}".format(**self.foopaths).split()
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args)

    def test_arg_parser_taxassign(self):

        # Ok
        args = "taxassign --asvtable {asvtable_tsv} --output {foodir} --blastdbdir {foodir} --blastdbname {filenottsv} " \
               "--taxonomy  {taxonomy_tsv}".format(**self.foopaths).split()
        self.assertTrue(self.parser.parse_args(args), 0)

        #######################################################################
        #
        # raises SystemExit
        #
        #######################################################################

        args = ["taxassign"]
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args)

    def test_arg_parser_taxassign_wrong_header_sequences(self):

        #######################################################################
        #
        # Wrong
        #
        #######################################################################

        self.foopaths['asvtable_tsv_wrong_header_sequences'] =  os.path.join(
            self.test_path, "test_files_dryad.f40v5_small", "run1_mfzr_zfzr", "asvtable_default_wrong_header_sequence.tsv")
        args = "taxassign --asvtable {asvtable_tsv_wrong_header_sequences} --output {foodir} --blastdbdir {foodir} --blastdbname {filenottsv} " \
               "--taxonomy  {taxonomy_tsv}".format(**self.foopaths).split()

        with self.assertRaises(SystemExit):
            self.assertTrue(self.parser.parse_args(args), 0)

    def test_arg_parser_pool_runmarker(self):

        # Ok
        args = "pool --db {filenottsv} --runmarker {runmarker_tsv} --asvtable {asvtable_tsv}".format(
            **self.foopaths).split()
        self.assertTrue(self.parser.parse_args(args), 0)

        #######################################################################
        #
        # raises SystemExit
        #
        #######################################################################

        args = ["pool"]
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args)

    def tearDown(self):
        shutil.rmtree(self.foopaths['sorteddir'], ignore_errors=True)
