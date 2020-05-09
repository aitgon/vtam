# -*- coding: utf-8 -*-
import os
import shutil

from pathlib import Path
from unittest import TestCase
from vtam.utils.ArgParser import ArgParser
from vtam.utils.PathManager import PathManager


class TestArgParser(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.parser = ArgParser.get_main_arg_parser()

        foopaths = {}
        foopaths['filedoesnotexist'] = "filedoesnotexist"
        foopaths['dirdoesnotexist'] = "dirdoesnotexist"
        foopaths['fileisempty'] = os.path.relpath("../test_files/emptyfile", PathManager.get_package_path())
        foopaths['filenottsv'] = os.path.relpath(__file__, PathManager.get_package_path())
        foopaths['readinfo_tsv'] = os.path.relpath(os.path.join(PathManager.get_test_path(), "test_files",
                                                                "readinfo.tsv"), PathManager.get_package_path())
        foopaths['params_yml'] = os.path.relpath(os.path.join(PathManager.get_test_path(), "test_files",
                                                                "params.yml"), PathManager.get_package_path())
        foopaths['params_wrong_yml'] = os.path.relpath(os.path.join(PathManager.get_test_path(), "test_files",
                                                                "params_wrong.yml"), PathManager.get_package_path())
        foopaths['fastqinfo_tsv_path'] = os.path.relpath(os.path.join(PathManager.get_test_path(), "test_files",
                                                                "known_occurrences.tsv"), PathManager.get_package_path())
        foopaths['asvtable_tsv'] = os.path.relpath(os.path.join(PathManager.get_test_path(), "test_files",
                                                                "asvtable.tsv"), PathManager.get_package_path())
        foopaths['runmarker_tsv'] = os.path.relpath(os.path.join(PathManager.get_test_path(), "test_files",
                                                                "pool_run_marker.tsv"), PathManager.get_package_path())
        foopaths['taxonomy_tsv'] = os.path.relpath(os.path.join(PathManager.get_test_path(), "test_files",
                                                                "taxonomy.tsv"), PathManager.get_package_path())
        foopaths['foodir'] = os.path.relpath(os.path.dirname(__file__), PathManager.get_package_path())
        foopaths['outdir'] = os.path.relpath(os.path.join(PathManager.get_test_path(),
                                                                             'output'), PathManager.get_package_path())
        foopaths['emptydir'] = os.path.relpath(os.path.join(foopaths['outdir'], 'emptydir'),
                                               PathManager.get_package_path())
        Path(os.path.join(foopaths['emptydir'])).mkdir(parents=True, exist_ok=True)
        foopaths['blastdb'] = os.path.relpath(os.path.join(PathManager.get_test_path(), 'test_files', 'blastdb'),
                                              PathManager.get_package_path())
        foopaths['known_occurrences_tsv'] = os.path.relpath(os.path.join(
            PathManager.get_package_path(), 'doc/data/dryad.f40v5_small/known_occurrences.tsv'),
                                              PathManager.get_package_path())
        cls.foopaths = foopaths

    def test_arg_parser_params(self):

        args = "filter --readinfo {readinfo_tsv} --readdir {foodir} --asvtable asvtable.tsv --params {params_yml}".format(
            **self.foopaths).split()
        self.assertTrue(self.parser.parse_args(args), 0)
        # Eror
        # import pdb; pdb.set_trace()
        args = "filter --readinfo {readinfo_tsv} --readdir {foodir} --asvtable asvtable.tsv --params {params_wrong_yml}".format(
            **self.foopaths).split()
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args)

    def test_arg_parser_filter(self):

        # Ok
        # import pdb; pdb.set_trace()
        args = "filter --readinfo {readinfo_tsv} --readdir {foodir} --asvtable asvtable.tsv".format(**self.foopaths).split()
        self.assertTrue(self.parser.parse_args(args), 0)

        ################################################################################################################
        #
        # raises SystemExit
        #
        ################################################################################################################

        args = ["filter"]
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args)

        args = "filter --readinfo {filedoesnotexist} --readdir {foodir} --outdir {foodir}".format(**self.foopaths).split()
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args)

        args = "filter --readinfo {filenottsv} --readdir {foodir} --outdir {foodir}".format(**self.foopaths).split()
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args)

        args = "filter --readinfo {readinfo_tsv} --readdir {emptydir} --outdir {foodir}".format(**self.foopaths).split()
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args)

        args = "filter --readinfo {readinfo_tsv} --readdir {dirdoesnotexist} --outdir {foodir}".format(**self.foopaths).split()
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args)

    def test_arg_parser_optimize(self):

        # Ok
        # TODO fix Mai 9, 2020
        args = "optimize --known_occurrences {known_occurrences_tsv} --readinfo {readinfo_tsv} --readdir {foodir} " \
               "--outdir {foodir}".format(**self.foopaths).split()
        self.assertTrue(self.parser.parse_args(args), 0)

        ################################################################################################################
        #
        # raises SystemExit
        #
        ################################################################################################################

        args = ["optimize"]
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args)

        args = "optimize --known_occurrences {filedoesnotexist} --readinfo {readinfo_tsv} --readdir {foodir} " \
               "--outdir {foodir}".format(**self.foopaths).split()
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args)

        args = "optimize --known_occurrences {filenottsv} --readinfo {readinfo_tsv} --readdir {foodir} " \
               "--outdir {foodir}".format(**self.foopaths).split()
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args)

    def test_arg_parser_taxassign(self):

        # Ok
        args = "taxassign --variants {asvtable_tsv} --output {foodir} --blastdbdir {foodir} --blastdbname {filenottsv} " \
               "--taxonomy  {taxonomy_tsv}".format(**self.foopaths).split()
        self.assertTrue(self.parser.parse_args(args), 0)

        ################################################################################################################
        #
        # raises SystemExit
        #
        ################################################################################################################

        args = ["taxassign"]
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args)

    def test_arg_parser_pool_runmarker(self):

        # Ok
        args = "pool --db {filenottsv} --runmarker {runmarker_tsv} --output {outdir}".format(**self.foopaths).split()
        self.assertTrue(self.parser.parse_args(args), 0)

        ################################################################################################################
        #
        # raises SystemExit
        #
        ################################################################################################################

        args = ["pool"]
        with self.assertRaises(SystemExit):
            self.parser.parse_args(args)

    def tearDown(self):
        shutil.rmtree(self.foopaths['outdir'], ignore_errors=True)
