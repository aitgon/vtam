import inspect
import os
import pathlib
from unittest import TestCase

import yaml

from vtam.utils.ArgParser import ArgParser

from vtam.utils.OptionManager import OptionManager
from vtam.utils.WopmarsRunner import WopmarsRunner
from vtam.utils.PathManager import PathManager


class TestWorpmarsRunnerTaxassign(TestCase):

    @classmethod
    def setUpClass(cls):

        foopaths = {}
        foopaths['foofile'] = os.path.relpath(__file__, PathManager.get_package_path())
        foopaths['foodir'] = os.path.relpath(os.path.dirname(__file__), PathManager.get_package_path())
        foopaths['outdir'] = os.path.relpath(os.path.join(PathManager.get_module_test_path(),
                                                                             'output'), PathManager.get_package_path())
        foopaths['blastdb'] = os.path.relpath(os.path.join(PathManager.get_module_test_path(), 'test_files', 'blastdb'),
                                              PathManager.get_package_path())
        cls.foopaths = foopaths

        cls.minseqlength_value_32 = 32
        cls.minseqlength_value_40 = 40
        cls.lfn_variant_replicate_threshold = 0.002

    def setUp(self):
        OptionManager.instance().clear()
        self.tempdir = PathManager.instance().get_tempdir()
        pathlib.Path(self.tempdir).mkdir(parents=True, exist_ok=True)

    def test_wopmars_runner_taxassign(self):
        #
        args_str = 'taxassign --taxonomy {foofile} --blast_db {blastdb}'.format(**self.foopaths)
        parser = ArgParser.get_arg_parser(is_abspath=False)
        args = parser.parse_args(args_str.split())

        #####################
        #
        # Add argparser attributes to optionmanager
        #
        #####################

        OptionManager.instance().clear()
        option_dic = vars(args)  # Dictionnary with options
        OptionManager.instance().add_options(option_dic)  # Add options to OptionManager

        ###############################################################
        #
        # Test wopfile
        #
        ###############################################################
        wopmars_runner = WopmarsRunner(command='taxassign', parameters=OptionManager.instance())
        wopfile_path, wopfile_content = wopmars_runner.create_wopfile()
        wopfile_content_bak = """rule TaxAssign:
    tool: vtam.wrapper.TaxAssign
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Variant: vtam.models.Variant
            FilterCodonStop: vtam.models.FilterCodonStop
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
            taxonomy: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            TaxAssign: vtam.models.TaxAssign
    params:
        update_taxassign: 0
        ltg_rule_threshold: 97
        include_prop: 90
        min_number_of_taxa: 3
        blast_db: test/test_files/blastdb"""
        import pdb; pdb.set_trace()
        self.assertTrue(wopfile_content == wopfile_content_bak)
