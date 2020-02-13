import os
from unittest import TestCase
from vtam.utils.ArgParser import ArgParser

from vtam.utils.OptionManager import OptionManager
from vtam.utils.WopmarsRunner import WopmarsRunner
from vtam.utils.PathManager import PathManager


class TestWorpmarsRunnerOptimize(TestCase):

    def setUp(self):
        # Minimal merge command
        foopaths = {}
        foopaths['foofile'] = os.path.relpath(__file__, PathManager.get_package_path())
        foopaths['foodir'] = os.path.relpath(os.path.dirname(__file__), PathManager.get_package_path())
        foopaths['outdir'] = os.path.relpath(os.path.join(PathManager.get_module_test_path(),
                                                                             'output'), PathManager.get_package_path())
        foopaths['blastdb'] = os.path.relpath(os.path.join(PathManager.get_module_test_path(), 'test_files', 'blastdb'),
                                              PathManager.get_package_path())
        self.foopaths = foopaths


    def test_wopmars_runner_optimize(self):
        args_str = 'optimize --fastainfo {foofile} --fastadir {outdir} --variant_known {foofile} --outdir {outdir}'\
            .format(**self.foopaths)
        parser = ArgParser.get_main_arg_parser()
        # import pdb; pdb.set_trace()
        args = parser.parse_args(args_str.split())

        #####################
        #
        # Add argparser attributes to optionmanager
        #
        #####################

        option_dic = vars(args) # Dictionnary with options
        OptionManager.instance().add_options(option_dic) # Add options to OptionManager

        ###############################################################
        #
        # Test wopfile
        #
        ###############################################################
        wopmars_runner = WopmarsRunner(command='optimize', parameters=OptionManager.instance())
        wopfile_path = os.path.relpath(os.path.join(PathManager.get_package_path(), "test/output/wopfile"),
                                    PathManager.get_package_path())
        wopfile_path, wopfile_content = wopmars_runner.create_wopfile(path=wopfile_path)
        wopfile_content_bak = """rule SampleInformation:
    tool: vtam.wrapper.SampleInformation
    input:
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_optimize.py
    output:
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
            Fasta: vtam.models.Fasta
            PrimerPair: vtam.models.PrimerPair
            TagPair: vtam.models.TagPair
            SampleInformation: vtam.models.SampleInformation


rule SortReads:
    tool: vtam.wrapper.SortReads
    input:
        table:
            Fasta: vtam.models.Fasta
            SampleInformation: vtam.models.SampleInformation
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_optimize.py
    output:
        file:
            sortreads: test/output/sortreads.tsv
    params:
        min_id: 0.8
        minseqlength: 32
        overhang: 0


rule VariantReadCount:
    tool: vtam.wrapper.VariantReadCount
    input:
        file:
            sortreads: test/output/sortreads.tsv
            fastainfo: test/utils/test_wopmars_runner_wopfile_optimize.py
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
    output:
        table:
            Variant: vtam.models.Variant
            VariantReadCount: vtam.models.VariantReadCount


rule OptimizeLFNbiosampleReplicate:
    tool: vtam.wrapper.OptimizeLFNbiosampleReplicate
    input:
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
            Variant: vtam.models.Variant
            VariantReadCount: vtam.models.VariantReadCount
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_optimize.py
            variant_known: test/utils/test_wopmars_runner_wopfile_optimize.py
    output:
        file:
            optimize_lfn_biosample_replicate: test/output/optimize_lfn_biosample_replicate.tsv


rule OptimizePCRerror:
    tool: vtam.wrapper.OptimizePCRerror
    input:
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
            Variant: vtam.models.Variant
            VariantReadCount: vtam.models.VariantReadCount
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_optimize.py
            variant_known: test/utils/test_wopmars_runner_wopfile_optimize.py
    output:
        file:
            optimize_pcr_error: test/output/optimize_pcr_error.tsv


rule OptimizeLFNreadCountAndLFNvariant:
    tool: vtam.wrapper.OptimizeLFNreadCountAndLFNvariant
    input:
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
            Variant: vtam.models.Variant
            VariantReadCount: vtam.models.VariantReadCount
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_optimize.py
            variant_known: test/utils/test_wopmars_runner_wopfile_optimize.py
    output:
        file:
            optimize_lfn_read_count_and_lfn_variant: test/output/optimize_lfn_read_count_and_lfn_variant.tsv
            optimize_lfn_variant_specific: test/output/optimize_lfn_variant_specific.tsv
    params:
        is_optimize_lfn_variant_replicate: 0
        lfn_variant_or_variant_replicate_threshold: 0.001
        lfn_biosample_replicate_threshold: 0.001
        lfn_read_count_threshold: 10
        min_replicate_number: 2"""
        self.assertTrue(wopfile_content == wopfile_content_bak)
