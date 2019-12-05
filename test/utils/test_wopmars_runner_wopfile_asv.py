import os
from unittest import TestCase
from vtam.utils.ArgParser import ArgParser

from vtam.utils.OptionManager import OptionManager
from vtam.utils.WopmarsRunner import WopmarsRunner
from vtam.utils.PathManager import PathManager


class TestWorpmarsRunnerASV(TestCase):

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


    def test_wopmars_runner_asv(self):
        #
        args_str = 'asv --fastainfo {foofile} --fastadir {foodir} --outdir {outdir} --taxonomy {foofile} --blast_db ' \
                   '{blastdb}'.format(**self.foopaths)
        parser = ArgParser.get_arg_parser(is_abspath=False)
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
        wopmars_runner = WopmarsRunner(command='asv', parameters=OptionManager.instance())
        wopfile_path, wopfile_content = wopmars_runner.create_wopfile()
        wopfile_content_bak = """rule SampleInformation:
    tool: vtam.wrapper.SampleInformation
    input:
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
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
            Replicate: vtam.models.Replicate
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
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
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
    output:
        table:
            Variant: vtam.models.Variant
            VariantReadCount: vtam.models.VariantReadCount


rule FilterLFN:
    tool: vtam.wrapper.FilterLFN
    input:
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            VariantReadCount: vtam.models.VariantReadCount
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            FilterLFN: vtam.models.FilterLFN
    params:
        filter_lfn_variant: 1
        lfn_variant_threshold: 0.001
        lfn_variant_replicate_threshold: 0.001
        lfn_biosample_replicate_threshold: 0.001
        lfn_read_count_threshold: 10


rule FilterMinReplicateNumber:
    tool: vtam.wrapper.FilterMinReplicateNumber
    input:
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            FilterLFN: vtam.models.FilterLFN
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
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
            Replicate: vtam.models.Replicate
            Variant: vtam.models.Variant
            FilterMinReplicateNumber: vtam.models.FilterMinReplicateNumber
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
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
            Replicate: vtam.models.Replicate
            Variant: vtam.models.Variant
            FilterPCRerror: vtam.models.FilterPCRerror
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            FilterChimera: vtam.models.FilterChimera
            FilterChimeraBorderline: vtam.models.FilterChimeraBorderline


rule FilterRenkonen:
    tool: vtam.wrapper.FilterRenkonen
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            FilterChimera: vtam.models.FilterChimera
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            FilterRenkonen: vtam.models.FilterRenkonen
    params:
        upper_renkonen_tail: 0.1


rule FilterIndel:
    tool: vtam.wrapper.FilterIndel
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            Variant: vtam.models.Variant
            FilterRenkonen: vtam.models.FilterRenkonen
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            FilterIndel: vtam.models.FilterIndel


rule FilterCodonStop:
    tool: vtam.wrapper.FilterCodonStop
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            Variant: vtam.models.Variant
            FilterIndel: vtam.models.FilterIndel
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            FilterCodonStop: vtam.models.FilterCodonStop
    params:
        genetic_table_number: 5


rule ReadCountAverageOverReplicates:
    tool: vtam.wrapper.ReadCountAverageOverReplicates
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            FilterCodonStop: vtam.models.FilterCodonStop
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            ReadCountAverageOverReplicates: vtam.models.ReadCountAverageOverReplicates


rule TaxAssign:
    tool: vtam.wrapper.TaxAssign
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            Variant: vtam.models.Variant
            FilterCodonStop: vtam.models.FilterCodonStop
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
            taxonomy: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            TaxAssign: vtam.models.TaxAssign
    params:
        ltg_rule_threshold: 97
        include_prop: 90
        min_number_of_taxa: 3
        blast_db: test/test_files/blastdb


rule MakeAsvTable:
    tool: vtam.wrapper.MakeAsvTable
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            Variant: vtam.models.Variant
            FilterChimeraBorderline: vtam.models.FilterChimeraBorderline
            FilterCodonStop: vtam.models.FilterCodonStop
            TaxAssign: vtam.models.TaxAssign
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
            taxonomy: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        file:
            ASVTable: test/output/asvtable.tsv


rule PoolMarkers:
    tool: vtam.wrapper.PoolMarkers
    input:
        file:
            ASVtable: test/output/asvtable.tsv
    output:
        file:
            PooledMarkers: test/output/pooled_markers.tsv"""
        self.assertTrue(wopfile_content == wopfile_content_bak)

    def test_wopmars_runner_asv_with_map_taxids(self):
            #
            # Minimal merge command
            # foopaths = {}
            # foopaths['foofile'] = os.path.relpath(__file__, PathManager.get_package_path())
            # foopaths['foodir'] = os.path.relpath(os.path.dirname(__file__), PathManager.get_package_path())
            # foopaths['outdir'] = os.path.join(PathManager.get_module_test_path(), 'output')
            # foopaths['blastdb'] = os.path.relpath(
            #     os.path.join(PathManager.get_module_test_path(), 'test_files', 'blastdb'),
            #     PathManager.get_package_path())
            # args_str = 'asv --fastqinfo {foofile} --fastqdir {foodir} --fastainfo {foofile} --fastadir {foodir}'.format(**foopaths)
            args_str = 'asv --fastainfo {foofile} --fastadir {foodir} --outdir {outdir} --taxonomy {foofile} --blast_db ' \
                       '{blastdb}'.format(**self.foopaths)
            parser = ArgParser.get_arg_parser(is_abspath=False)
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
            wopmars_runner = WopmarsRunner(command='asv', parameters=OptionManager.instance())
            wopfile_path = os.path.relpath(os.path.join(PathManager.get_package_path(), "test/output/wopfile"),
                                           PathManager.get_package_path())
            wopfile_path, wopfile_content = wopmars_runner.create_wopfile(path=wopfile_path)
            wopfile_content_bak = """rule SampleInformation:
    tool: vtam.wrapper.SampleInformation
    input:
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
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
            Replicate: vtam.models.Replicate
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
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
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
    output:
        table:
            Variant: vtam.models.Variant
            VariantReadCount: vtam.models.VariantReadCount


rule FilterLFN:
    tool: vtam.wrapper.FilterLFN
    input:
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            VariantReadCount: vtam.models.VariantReadCount
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            FilterLFN: vtam.models.FilterLFN
    params:
        filter_lfn_variant: 1
        lfn_variant_threshold: 0.001
        lfn_variant_replicate_threshold: 0.001
        lfn_biosample_replicate_threshold: 0.001
        lfn_read_count_threshold: 10


rule FilterMinReplicateNumber:
    tool: vtam.wrapper.FilterMinReplicateNumber
    input:
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            FilterLFN: vtam.models.FilterLFN
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
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
            Replicate: vtam.models.Replicate
            Variant: vtam.models.Variant
            FilterMinReplicateNumber: vtam.models.FilterMinReplicateNumber
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
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
            Replicate: vtam.models.Replicate
            Variant: vtam.models.Variant
            FilterPCRerror: vtam.models.FilterPCRerror
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            FilterChimera: vtam.models.FilterChimera
            FilterChimeraBorderline: vtam.models.FilterChimeraBorderline


rule FilterRenkonen:
    tool: vtam.wrapper.FilterRenkonen
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            FilterChimera: vtam.models.FilterChimera
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            FilterRenkonen: vtam.models.FilterRenkonen
    params:
        upper_renkonen_tail: 0.1


rule FilterIndel:
    tool: vtam.wrapper.FilterIndel
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            Variant: vtam.models.Variant
            FilterRenkonen: vtam.models.FilterRenkonen
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            FilterIndel: vtam.models.FilterIndel


rule FilterCodonStop:
    tool: vtam.wrapper.FilterCodonStop
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            Variant: vtam.models.Variant
            FilterIndel: vtam.models.FilterIndel
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            FilterCodonStop: vtam.models.FilterCodonStop
    params:
        genetic_table_number: 5


rule ReadCountAverageOverReplicates:
    tool: vtam.wrapper.ReadCountAverageOverReplicates
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            FilterCodonStop: vtam.models.FilterCodonStop
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            ReadCountAverageOverReplicates: vtam.models.ReadCountAverageOverReplicates


rule TaxAssign:
    tool: vtam.wrapper.TaxAssign
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            Variant: vtam.models.Variant
            FilterCodonStop: vtam.models.FilterCodonStop
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
            taxonomy: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            TaxAssign: vtam.models.TaxAssign
    params:
        ltg_rule_threshold: 97
        include_prop: 90
        min_number_of_taxa: 3
        blast_db: test/test_files/blastdb


rule MakeAsvTable:
    tool: vtam.wrapper.MakeAsvTable
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            Variant: vtam.models.Variant
            FilterChimeraBorderline: vtam.models.FilterChimeraBorderline
            FilterCodonStop: vtam.models.FilterCodonStop
            TaxAssign: vtam.models.TaxAssign
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
            taxonomy: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        file:
            ASVTable: test/output/asvtable.tsv


rule PoolMarkers:
    tool: vtam.wrapper.PoolMarkers
    input:
        file:
            ASVtable: test/output/asvtable.tsv
    output:
        file:
            PooledMarkers: test/output/pooled_markers.tsv"""
            self.assertTrue(wopfile_content == wopfile_content_bak)




    def test_wopmars_runner_asv_with_threshold_specific(self):

        args_str = 'asv --fastainfo {foofile} --fastadir {foodir} --outdir {outdir} --taxonomy {foofile} --blast_db ' \
                   '{blastdb} --threshold_specific {foofile}'.format(**self.foopaths)
        parser = ArgParser.get_arg_parser(is_abspath=False)
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
        wopmars_runner = WopmarsRunner(command='asv', parameters=OptionManager.instance())
        wopfile_path = os.path.relpath(os.path.join(PathManager.get_package_path(), "test/output/wopfile"),
                                    PathManager.get_package_path())
        wopfile_path, wopfile_content = wopmars_runner.create_wopfile(path=wopfile_path)
        wopfile_content_bak = """rule SampleInformation:
    tool: vtam.wrapper.SampleInformation
    input:
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
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
            Replicate: vtam.models.Replicate
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
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
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
    output:
        table:
            Variant: vtam.models.Variant
            VariantReadCount: vtam.models.VariantReadCount


rule FilterLFN:
    tool: vtam.wrapper.FilterLFN
    input:
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            VariantReadCount: vtam.models.VariantReadCount
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
            threshold_specific: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            FilterLFN: vtam.models.FilterLFN
    params:
        filter_lfn_variant: 1
        lfn_variant_threshold: 0.001
        lfn_variant_replicate_threshold: 0.001
        lfn_biosample_replicate_threshold: 0.001
        lfn_read_count_threshold: 10


rule FilterMinReplicateNumber:
    tool: vtam.wrapper.FilterMinReplicateNumber
    input:
        table:
            Run: vtam.models.Run
            Marker: vtam.models.Marker
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            FilterLFN: vtam.models.FilterLFN
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
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
            Replicate: vtam.models.Replicate
            Variant: vtam.models.Variant
            FilterMinReplicateNumber: vtam.models.FilterMinReplicateNumber
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
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
            Replicate: vtam.models.Replicate
            Variant: vtam.models.Variant
            FilterPCRerror: vtam.models.FilterPCRerror
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            FilterChimera: vtam.models.FilterChimera
            FilterChimeraBorderline: vtam.models.FilterChimeraBorderline


rule FilterRenkonen:
    tool: vtam.wrapper.FilterRenkonen
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            FilterChimera: vtam.models.FilterChimera
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            FilterRenkonen: vtam.models.FilterRenkonen
    params:
        upper_renkonen_tail: 0.1


rule FilterIndel:
    tool: vtam.wrapper.FilterIndel
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            Variant: vtam.models.Variant
            FilterRenkonen: vtam.models.FilterRenkonen
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            FilterIndel: vtam.models.FilterIndel


rule FilterCodonStop:
    tool: vtam.wrapper.FilterCodonStop
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            Variant: vtam.models.Variant
            FilterIndel: vtam.models.FilterIndel
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            FilterCodonStop: vtam.models.FilterCodonStop
    params:
        genetic_table_number: 5


rule ReadCountAverageOverReplicates:
    tool: vtam.wrapper.ReadCountAverageOverReplicates
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            FilterCodonStop: vtam.models.FilterCodonStop
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            ReadCountAverageOverReplicates: vtam.models.ReadCountAverageOverReplicates


rule TaxAssign:
    tool: vtam.wrapper.TaxAssign
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            Variant: vtam.models.Variant
            FilterCodonStop: vtam.models.FilterCodonStop
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
            taxonomy: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        table:
            TaxAssign: vtam.models.TaxAssign
    params:
        ltg_rule_threshold: 97
        include_prop: 90
        min_number_of_taxa: 3
        blast_db: test/test_files/blastdb


rule MakeAsvTable:
    tool: vtam.wrapper.MakeAsvTable
    input:
        table:
            Marker: vtam.models.Marker
            Run: vtam.models.Run
            Biosample: vtam.models.Biosample
            Replicate: vtam.models.Replicate
            Variant: vtam.models.Variant
            FilterChimeraBorderline: vtam.models.FilterChimeraBorderline
            FilterCodonStop: vtam.models.FilterCodonStop
            TaxAssign: vtam.models.TaxAssign
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_asv.py
            taxonomy: test/utils/test_wopmars_runner_wopfile_asv.py
    output:
        file:
            ASVTable: test/output/asvtable.tsv


rule PoolMarkers:
    tool: vtam.wrapper.PoolMarkers
    input:
        file:
            ASVtable: test/output/asvtable.tsv
    output:
        file:
            PooledMarkers: test/output/pooled_markers.tsv"""
        self.assertTrue(wopfile_content == wopfile_content_bak)
