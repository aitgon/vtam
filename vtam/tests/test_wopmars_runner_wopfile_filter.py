from unittest import TestCase
from vtam.utils.ArgParser import ArgParser
from vtam.utils.OptionManager import OptionManager
from vtam.utils.PathManager import PathManager
from vtam.utils.WopmarsRunner import WopmarsRunner
import os
import pathlib
import yaml


class TestWorpmarsRunnerFilter(TestCase):

    @classmethod
    def setUpClass(cls):

        foopaths = {}
        foopaths['foofile'] = os.path.relpath(__file__, PathManager.get_package_path())
        foopaths['foodir'] = os.path.relpath(os.path.dirname(__file__), PathManager.get_package_path())
        foopaths['outdir'] = os.path.relpath(os.path.join(PathManager.get_test_path(),
                                                                             'output'), PathManager.get_package_path())
        foopaths['blastdb'] = os.path.relpath(os.path.join(PathManager.get_test_path(), 'test_files', 'blastdb'),
                                              PathManager.get_package_path())
        foopaths['readinfo_tsv'] = os.path.relpath(os.path.join(PathManager.get_test_path(), "test_files",
                                                                "readinfo.tsv"), PathManager.get_package_path())
        cls.foopaths = foopaths

        cls.minseqlength_value_32 = 32
        cls.minseqlength_value_40 = 40
        cls.lfn_variant_replicate_threshold = 0.002

    def setUp(self):
        OptionManager.instance().clear()
        self.tempdir = PathManager.instance().get_tempdir()
        pathlib.Path(self.tempdir).mkdir(parents=True, exist_ok=True)

    def test_wopmars_runner_filter(self):
        #
        args_str = 'filter --readinfo {readinfo_tsv} --readdir {foodir} --outdir {outdir}'.format(**self.foopaths)
        parser = ArgParser.get_main_arg_parser()
        args = parser.parse_args(args_str.split())

        #####################
        #
        # Add argparser attributes to optionmanager
        #
        #####################

        OptionManager.instance().clear()
        option_dic = vars(args) # Dictionnary with options
        OptionManager.instance().add_options(option_dic) # Add options to OptionManager

        ###############################################################
        #
        # Test wopfile
        #
        ###############################################################
        wopmars_runner = WopmarsRunner(command='filter', parameters=OptionManager.instance())
        wopfile_path, wopfile_content = wopmars_runner.create_wopfile()
        wopfile_content_bak = """rule SampleInformation:
    tool: vtam.wrapper.SampleInformation
    input:
        file:
            readinfo: vtam/tests/test_files/readinfo.tsv
    output:
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
            Fasta: vtam.models.Fasta
            SampleInformation: vtam.models.SampleInformation


rule VariantReadCount:
    tool: vtam.wrapper.VariantReadCount
    input:
        file:
            readinfo: vtam/tests/test_files/readinfo.tsv
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
    output:
        table:
            Variant: vtam.models.Variant
            VariantReadCount: vtam.models.VariantReadCount
    params:
        read_dir: vtam/tests


rule FilterLFN:
    tool: vtam.wrapper.FilterLFN
    input:
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
            VariantReadCount: vtam.models.VariantReadCount
        file:
            readinfo: vtam/tests/test_files/readinfo.tsv
    output:
        table:
            FilterLFN: vtam.models.FilterLFN
    params:
        lfn_variant_threshold: 0.001
        lfn_biosample_replicate_threshold: 0.001
        lfn_read_count_threshold: 10


rule FilterMinReplicateNumber:
    tool: vtam.wrapper.FilterMinReplicateNumber
    input:
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
            FilterLFN: vtam.models.FilterLFN
        file:
            readinfo: vtam/tests/test_files/readinfo.tsv
    output:
        table:
            FilterMinReplicateNumber: vtam.models.FilterMinReplicateNumber
    params:
        min_replicate_number: 2


rule FilterPCRerror:
    tool: vtam.wrapper.FilterPCRerror
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Variant: vtam.models.Variant
            FilterMinReplicateNumber: vtam.models.FilterMinReplicateNumber
        file:
            readinfo: vtam/tests/test_files/readinfo.tsv
    output:
        table:
            FilterPCRerror: vtam.models.FilterPCRerror
    params:
        pcr_error_var_prop: 0.1


rule FilterChimera:
    tool: vtam.wrapper.FilterChimera
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Variant: vtam.models.Variant
            FilterPCRerror: vtam.models.FilterPCRerror
        file:
            readinfo: vtam/tests/test_files/readinfo.tsv
    output:
        table:
            FilterChimera: vtam.models.FilterChimera
            FilterChimeraBorderline: vtam.models.FilterChimeraBorderline


rule FilterMinReplicateNumber2:
    tool: vtam.wrapper.FilterMinReplicateNumber
    input:
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
            FilterLFN: vtam.models.FilterChimera
        file:
            readinfo: vtam/tests/test_files/readinfo.tsv
    output:
        table:
            FilterMinReplicateNumber: vtam.models.FilterMinReplicateNumber2
    params:
        min_replicate_number: 2


rule FilterRenkonen:
    tool: vtam.wrapper.FilterRenkonen
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            FilterChimera: vtam.models.FilterMinReplicateNumber2
        file:
            readinfo: vtam/tests/test_files/readinfo.tsv
    output:
        table:
            FilterRenkonen: vtam.models.FilterRenkonen
    params:
        upper_renkonen_tail: 0.1


rule FilterMinReplicateNumber3:
    tool: vtam.wrapper.FilterMinReplicateNumber
    input:
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
            FilterLFN: vtam.models.FilterRenkonen
        file:
            readinfo: vtam/tests/test_files/readinfo.tsv
    output:
        table:
            FilterMinReplicateNumber: vtam.models.FilterMinReplicateNumber3
    params:
        min_replicate_number: 2


rule FilterIndel:
    tool: vtam.wrapper.FilterIndel
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Variant: vtam.models.Variant
            FilterRenkonen: vtam.models.FilterMinReplicateNumber3
        file:
            readinfo: vtam/tests/test_files/readinfo.tsv
    output:
        table:
            FilterIndel: vtam.models.FilterIndel
    params:
        skip_filter_indel: 0


rule FilterCodonStop:
    tool: vtam.wrapper.FilterCodonStop
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Variant: vtam.models.Variant
            FilterIndel: vtam.models.FilterIndel
        file:
            readinfo: vtam/tests/test_files/readinfo.tsv
    output:
        table:
            FilterCodonStop: vtam.models.FilterCodonStop
    params:
        genetic_table_number: 5
        skip_filter_codon_stop: 0


rule ReadCountAverageOverReplicates:
    tool: vtam.wrapper.ReadCountAverageOverReplicates
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            FilterCodonStop: vtam.models.FilterCodonStop
        file:
            readinfo: vtam/tests/test_files/readinfo.tsv
    output:
        table:
            ReadCountAverageOverReplicates: vtam.models.ReadCountAverageOverReplicates


rule MakeAsvTable:
    tool: vtam.wrapper.MakeAsvTable
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Variant: vtam.models.Variant
            FilterChimeraBorderline: vtam.models.FilterChimeraBorderline
            FilterCodonStop: vtam.models.FilterCodonStop
        file:
            readinfo: vtam/tests/test_files/readinfo.tsv
    output:
        file:
            ASVTable: vtam/tests/output/asvtable.tsv"""
        self.assertTrue(wopfile_content == wopfile_content_bak)

    def test_wopmars_runner_asv_with_threshold_specific(self):

        args_str = 'filter --readinfo {readinfo_tsv} --readdir {foodir} --outdir {outdir}' \
                   ' --threshold_specific {foofile}'.format(**self.foopaths)
        parser = ArgParser.get_main_arg_parser()
        args = parser.parse_args(args_str.split())

        #####################
        #
        # Add argparser attributes to optionmanager
        #
        #####################

        option_dic = vars(args)  # Dictionnary with options
        OptionManager.instance().add_options(option_dic)  # Add options to OptionManager

        ###############################################################
        #
        # Test wopfile
        #
        ###############################################################

        wopmars_runner = WopmarsRunner(command='filter', parameters=OptionManager.instance())
        wopfile_path = os.path.relpath(os.path.join(PathManager.get_package_path(), "tests/output/wopfile"),
                                    PathManager.get_package_path())
        wopfile_path, wopfile_content = wopmars_runner.create_wopfile(path=wopfile_path)

        self.assertTrue(yaml.load(wopfile_content, Loader=yaml.SafeLoader)['rule FilterLFN']['input']['file']['threshold_specific']
                        == self.foopaths['foofile'])

    def test_wopmars_runner_asv_with_lfn_variant_replicate(self):

        #####################
        #
        # Params yml
        #
        #####################

        params_yml_str = "lfn_variant_replicate_threshold: {}".format(self.lfn_variant_replicate_threshold)
        params_yml_path = os.path.join(self.tempdir, "params.yml")
        with open(params_yml_path, "w") as fout:
            fout.write(params_yml_str)
        this_foopaths = self.foopaths.copy()
        this_foopaths['params_yml'] = params_yml_path

        args_str = 'filter --readinfo {readinfo_tsv} --readdir {foodir} --outdir {outdir} --threshold_specific {foofile} ' \
                   '--params {params_yml}'.format(**this_foopaths)
        parser = ArgParser.get_main_arg_parser()
        args = parser.parse_args(args_str.split())

        #####################
        #
        # Add argparser attributes to optionmanager
        #
        #####################

        option_dic = vars(args)  # Dictionnary with options
        OptionManager.instance().add_options(option_dic)  # Add options to OptionManager

        ###############################################################
        #
        # Test wopfile
        #
        ###############################################################

        wopmars_runner = WopmarsRunner(command='filter', parameters=OptionManager.instance())
        wopfile_path = os.path.relpath(os.path.join(PathManager.get_package_path(), "tests/output/wopfile"),
                                    PathManager.get_package_path())
        wopfile_path, wopfile_content = wopmars_runner.create_wopfile(path=wopfile_path)

        self.assertTrue('lfn_variant_replicate_threshold' in yaml.load(wopfile_content, Loader=yaml.SafeLoader)['rule FilterLFN']['params'])
        self.assertFalse('lfn_variant_threshold' in yaml.load(wopfile_content, Loader=yaml.SafeLoader)['rule FilterLFN']['params'])
        self.assertTrue(yaml.load(wopfile_content, Loader=yaml.SafeLoader)['rule FilterLFN']['params']['lfn_variant_replicate_threshold']
                        == self.lfn_variant_replicate_threshold)