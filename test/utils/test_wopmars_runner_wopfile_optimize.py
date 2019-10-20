import os
from unittest import TestCase
from vtam.utils.ArgParser import ArgParser

from vtam.utils.OptionManager import OptionManager
from vtam.utils.WopmarsRunner import WopmarsRunner
from vtam.utils.PathManager import PathManager


class TestWorpmarsRunnerOptimize(TestCase):


    def test_wopmars_runner_optimize(self):
        #
        # Minimal merge command
        foopaths = {}
        foopaths['foofile'] = os.path.relpath(__file__, PathManager.get_package_path())
        foopaths['foodir'] = os.path.relpath(os.path.dirname(__file__), PathManager.get_package_path())
        args_str = 'optimize --fastainfo {foofile} --fastadir {foodir} --variant_known {foofile} --outdir {foodir}'.format(**foopaths)
        parser = ArgParser.get_arg_parser(is_abspath=False)
        # import pdb; pdb.set_trace()
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
            Run: vtam.model.Run
            Marker: vtam.model.Marker
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Fasta: vtam.model.Fasta
            PrimerPair: vtam.model.PrimerPair
            TagPair: vtam.model.TagPair
            SampleInformation: vtam.model.SampleInformation
    params:
        fasta_dir: test/utils


rule SortReads:
    tool: vtam.wrapper.SortReads
    input:
        table:
            Fasta: vtam.model.Fasta
            SampleInformation: vtam.model.SampleInformation
            Marker: vtam.model.Marker
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_optimize.py
    output:
        file:
            sortreads: test/utils/sortreads.tsv
    params:
        min_id: 0.8
        minseqlength: 32
        overhang: 0


rule VariantReadCount:
    tool: vtam.wrapper.VariantReadCount
    input:
        file:
            sortreads: test/utils/sortreads.tsv
            fastainfo: test/utils/test_wopmars_runner_wopfile_optimize.py
        table:
            Run: vtam.model.Run
            Marker: vtam.model.Marker
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
    output:
        table:
            Variant: vtam.model.Variant
            VariantReadCount: vtam.model.VariantReadCount
    params:
        foo: 0


rule OptimizeLFNbiosampleReplicate:
    tool: vtam.wrapper.OptimizeLFNbiosampleReplicate
    input:
        table:
            Run: vtam.model.Run
            Marker: vtam.model.Marker
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Variant: vtam.model.Variant
            VariantReadCount: vtam.model.VariantReadCount
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_optimize.py
            variant_known: test/utils/test_wopmars_runner_wopfile_optimize.py
    output:
        file:
            optimize_lfn_biosample_replicate: test/utils/optimize_lfn_biosample_replicate.tsv
    params:
      foo: 0


rule OptimizePCRerror:
    tool: vtam.wrapper.OptimizePCRerror
    input:
        table:
            Run: vtam.model.Run
            Marker: vtam.model.Marker
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Variant: vtam.model.Variant
            VariantReadCount: vtam.model.VariantReadCount
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_optimize.py
            variant_known: test/utils/test_wopmars_runner_wopfile_optimize.py
    output:
        file:
            optimize_pcr_error: test/utils/optimize_pcr_error.tsv
    params:
        foo: 0


rule OptimizeLFNreadCountAndLFNvariant:
    tool: vtam.wrapper.OptimizeLFNreadCountAndLFNvariant
    input:
        table:
            Run: vtam.model.Run
            Marker: vtam.model.Marker
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Variant: vtam.model.Variant
            VariantReadCount: vtam.model.VariantReadCount
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_optimize.py
            variant_known: test/utils/test_wopmars_runner_wopfile_optimize.py
    output:
        file:
            optimize_lfn_read_count_and_lfn_variant: test/utils/optimize_lfn_read_count_and_lfn_variant.tsv
            optimize_lfn_variant_specific: test/utils/optimize_lfn_variant_specific.tsv
    params:
        is_optimize_lfn_variant_replicate: 0
        lfn_variant_or_variant_replicate_threshold: 0.001
        lfn_biosample_replicate_threshold: 0.001
        lfn_read_count_threshold: 10
        min_replicate_number: 2


rule OptimizeLFNreadCountAndLFNvariantReplicate:
    tool: vtam.wrapper.OptimizeLFNreadCountAndLFNvariant
    input:
        table:
            Run: vtam.model.Run
            Marker: vtam.model.Marker
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Variant: vtam.model.Variant
            VariantReadCount: vtam.model.VariantReadCount
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_optimize.py
            variant_known: test/utils/test_wopmars_runner_wopfile_optimize.py
    output:
        file:
            optimize_lfn_read_count_and_lfn_variant: test/utils/optimize_lfn_read_count_and_lfn_variant_replicate.tsv
            optimize_lfn_variant_specific: test/utils/optimize_lfn_variant_replicate_specific.tsv
    params:
        is_optimize_lfn_variant_replicate: 1
        lfn_variant_or_variant_replicate_threshold: 0.001
        lfn_biosample_replicate_threshold: 0.001
        lfn_read_count_threshold: 10
        min_replicate_number: 2"""
        self.assertTrue(wopfile_content == wopfile_content_bak)
