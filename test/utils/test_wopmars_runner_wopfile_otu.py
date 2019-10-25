import os
from unittest import TestCase
from vtam.utils.ArgParser import ArgParser

from vtam.utils.OptionManager import OptionManager
from vtam.utils.WopmarsRunner import WopmarsRunner
from vtam.utils.PathManager import PathManager


class TestWorpmarsRunnerOTU(TestCase):

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


    def test_wopmars_runner_otu(self):
        #
        args_str = 'otu --fastainfo {foofile} --fastadir {foodir} --outdir {outdir} --taxonomy {foofile} --blast_db ' \
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
        wopmars_runner = WopmarsRunner(command='otu', parameters=OptionManager.instance())
        wopfile_path, wopfile_content = wopmars_runner.create_wopfile()
        wopfile_content_bak = """rule SampleInformation:
    tool: vtam.wrapper.SampleInformation
    input:
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
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
        log_verbosity: 0


rule SortReads:
    tool: vtam.wrapper.SortReads
    input:
        table:
            Fasta: vtam.model.Fasta
            SampleInformation: vtam.model.SampleInformation
            Marker: vtam.model.Marker
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        file:
            sortreads: test/output/sortreads.tsv
    params:
        min_id: 0.8
        minseqlength: 32
        overhang: 0
        log_verbosity: 0


rule VariantReadCount:
    tool: vtam.wrapper.VariantReadCount
    input:
        file:
            sortreads: test/output/sortreads.tsv
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
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
        log_verbosity: 0


rule FilterLFN:
    tool: vtam.wrapper.FilterLFN
    input:
        table:
            Run: vtam.model.Run
            Marker: vtam.model.Marker
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            VariantReadCount: vtam.model.VariantReadCount
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            FilterLFN: vtam.model.FilterLFN
    params:
        filter_lfn_variant: 1
        lfn_variant_threshold: 0.001
        lfn_variant_replicate_threshold: 0.001
        lfn_biosample_replicate_threshold: 0.001
        lfn_read_count_threshold: 10
        log_verbosity: 0


rule FilterMinReplicateNumber:
    tool: vtam.wrapper.FilterMinReplicateNumber
    input:
        table:
            Run: vtam.model.Run
            Marker: vtam.model.Marker
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            FilterLFN: vtam.model.FilterLFN
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            FilterMinReplicateNumber: vtam.model.FilterMinReplicateNumber
    params:
        min_replicate_number: 2
        log_verbosity: 0


rule FilterPCRError:
    tool: vtam.wrapper.FilterPCRError
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Variant: vtam.model.Variant
            FilterMinReplicateNumber: vtam.model.FilterMinReplicateNumber
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            FilterPCRError: vtam.model.FilterPCRError
    params:
        pcr_error_var_prop: 0.1
        log_verbosity: 0


rule FilterChimera:
    tool: vtam.wrapper.FilterChimera
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Variant: vtam.model.Variant
            FilterPCRError: vtam.model.FilterPCRError
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            FilterChimera: vtam.model.FilterChimera
            FilterChimeraBorderline: vtam.model.FilterChimeraBorderline
    params:
        log_verbosity: 0


rule FilterRenkonen:
    tool: vtam.wrapper.FilterRenkonen
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            FilterChimera: vtam.model.FilterChimera
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            FilterRenkonen: vtam.model.FilterRenkonen
    params:
        renkonen_threshold: 0.1
        log_verbosity: 0


rule FilterIndel:
    tool: vtam.wrapper.FilterIndel
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Variant: vtam.model.Variant
            FilterRenkonen: vtam.model.FilterRenkonen
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            FilterIndel: vtam.model.FilterIndel
    params:
        log_verbosity: 0


rule FilterCodonStop:
    tool: vtam.wrapper.FilterCodonStop
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Variant: vtam.model.Variant
            FilterIndel: vtam.model.FilterIndel
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            FilterCodonStop: vtam.model.FilterCodonStop
    params:
        genetic_table_number: 5
        log_verbosity: 0


rule ReadCountAverageOverReplicates:
    tool: vtam.wrapper.ReadCountAverageOverReplicates
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            FilterCodonStop: vtam.model.FilterCodonStop
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            ReadCountAverageOverReplicates: vtam.model.ReadCountAverageOverReplicates
    params:
        log_verbosity: 0


rule TaxAssign:
    tool: vtam.wrapper.TaxAssign
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Variant: vtam.model.Variant
            FilterCodonStop: vtam.model.FilterCodonStop
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
            taxonomy: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            TaxAssign: vtam.model.TaxAssign
    params:
        identity_threshold: 97
        include_prop: 90
        min_number_of_taxa: 3
        blast_db: test/test_files/blastdb
        num_threads: 8
        log_verbosity: 0


rule MakeOtuTable:
    tool: vtam.wrapper.MakeOtuTable
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Variant: vtam.model.Variant
            FilterChimeraBorderline: vtam.model.FilterChimeraBorderline
            FilterCodonStop: vtam.model.FilterCodonStop
            TaxAssign: vtam.model.TaxAssign
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
            taxonomy: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        file:
            OTUTable: test/output/otutable.tsv
    params:
        log_verbosity: 0"""
        self.assertTrue(wopfile_content == wopfile_content_bak)

    def test_wopmars_runner_otu_with_map_taxids(self):
            #
            # Minimal merge command
            # foopaths = {}
            # foopaths['foofile'] = os.path.relpath(__file__, PathManager.get_package_path())
            # foopaths['foodir'] = os.path.relpath(os.path.dirname(__file__), PathManager.get_package_path())
            # foopaths['outdir'] = os.path.join(PathManager.get_module_test_path(), 'output')
            # foopaths['blastdb'] = os.path.relpath(
            #     os.path.join(PathManager.get_module_test_path(), 'test_files', 'blastdb'),
            #     PathManager.get_package_path())
            # args_str = 'otu --fastqinfo {foofile} --fastqdir {foodir} --fastainfo {foofile} --fastadir {foodir}'.format(**foopaths)
            args_str = 'otu --fastainfo {foofile} --fastadir {foodir} --outdir {outdir} --taxonomy {foofile} --blast_db ' \
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
            wopmars_runner = WopmarsRunner(command='otu', parameters=OptionManager.instance())
            wopfile_path = os.path.relpath(os.path.join(PathManager.get_package_path(), "test/output/wopfile"),
                                           PathManager.get_package_path())
            wopfile_path, wopfile_content = wopmars_runner.create_wopfile(path=wopfile_path)
            wopfile_content_bak = """rule SampleInformation:
    tool: vtam.wrapper.SampleInformation
    input:
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
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
        log_verbosity: 0


rule SortReads:
    tool: vtam.wrapper.SortReads
    input:
        table:
            Fasta: vtam.model.Fasta
            SampleInformation: vtam.model.SampleInformation
            Marker: vtam.model.Marker
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        file:
            sortreads: test/output/sortreads.tsv
    params:
        min_id: 0.8
        minseqlength: 32
        overhang: 0
        log_verbosity: 0


rule VariantReadCount:
    tool: vtam.wrapper.VariantReadCount
    input:
        file:
            sortreads: test/output/sortreads.tsv
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
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
        log_verbosity: 0


rule FilterLFN:
    tool: vtam.wrapper.FilterLFN
    input:
        table:
            Run: vtam.model.Run
            Marker: vtam.model.Marker
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            VariantReadCount: vtam.model.VariantReadCount
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            FilterLFN: vtam.model.FilterLFN
    params:
        filter_lfn_variant: 1
        lfn_variant_threshold: 0.001
        lfn_variant_replicate_threshold: 0.001
        lfn_biosample_replicate_threshold: 0.001
        lfn_read_count_threshold: 10
        log_verbosity: 0


rule FilterMinReplicateNumber:
    tool: vtam.wrapper.FilterMinReplicateNumber
    input:
        table:
            Run: vtam.model.Run
            Marker: vtam.model.Marker
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            FilterLFN: vtam.model.FilterLFN
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            FilterMinReplicateNumber: vtam.model.FilterMinReplicateNumber
    params:
        min_replicate_number: 2
        log_verbosity: 0


rule FilterPCRError:
    tool: vtam.wrapper.FilterPCRError
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Variant: vtam.model.Variant
            FilterMinReplicateNumber: vtam.model.FilterMinReplicateNumber
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            FilterPCRError: vtam.model.FilterPCRError
    params:
        pcr_error_var_prop: 0.1
        log_verbosity: 0


rule FilterChimera:
    tool: vtam.wrapper.FilterChimera
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Variant: vtam.model.Variant
            FilterPCRError: vtam.model.FilterPCRError
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            FilterChimera: vtam.model.FilterChimera
            FilterChimeraBorderline: vtam.model.FilterChimeraBorderline
    params:
        log_verbosity: 0


rule FilterRenkonen:
    tool: vtam.wrapper.FilterRenkonen
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            FilterChimera: vtam.model.FilterChimera
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            FilterRenkonen: vtam.model.FilterRenkonen
    params:
        renkonen_threshold: 0.1
        log_verbosity: 0


rule FilterIndel:
    tool: vtam.wrapper.FilterIndel
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Variant: vtam.model.Variant
            FilterRenkonen: vtam.model.FilterRenkonen
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            FilterIndel: vtam.model.FilterIndel
    params:
        log_verbosity: 0


rule FilterCodonStop:
    tool: vtam.wrapper.FilterCodonStop
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Variant: vtam.model.Variant
            FilterIndel: vtam.model.FilterIndel
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            FilterCodonStop: vtam.model.FilterCodonStop
    params:
        genetic_table_number: 5
        log_verbosity: 0


rule ReadCountAverageOverReplicates:
    tool: vtam.wrapper.ReadCountAverageOverReplicates
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            FilterCodonStop: vtam.model.FilterCodonStop
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            ReadCountAverageOverReplicates: vtam.model.ReadCountAverageOverReplicates
    params:
        log_verbosity: 0


rule TaxAssign:
    tool: vtam.wrapper.TaxAssign
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Variant: vtam.model.Variant
            FilterCodonStop: vtam.model.FilterCodonStop
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
            taxonomy: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            TaxAssign: vtam.model.TaxAssign
    params:
        identity_threshold: 97
        include_prop: 90
        min_number_of_taxa: 3
        blast_db: test/test_files/blastdb
        num_threads: 8
        log_verbosity: 0


rule MakeOtuTable:
    tool: vtam.wrapper.MakeOtuTable
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Variant: vtam.model.Variant
            FilterChimeraBorderline: vtam.model.FilterChimeraBorderline
            FilterCodonStop: vtam.model.FilterCodonStop
            TaxAssign: vtam.model.TaxAssign
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
            taxonomy: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        file:
            OTUTable: test/output/otutable.tsv
    params:
        log_verbosity: 0"""
            self.assertTrue(wopfile_content == wopfile_content_bak)




    def test_wopmars_runner_otu_with_threshold_specific(self):

        args_str = 'otu --fastainfo {foofile} --fastadir {foodir} --outdir {outdir} --taxonomy {foofile} --blast_db ' \
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
        wopmars_runner = WopmarsRunner(command='otu', parameters=OptionManager.instance())
        wopfile_path = os.path.relpath(os.path.join(PathManager.get_package_path(), "test/output/wopfile"),
                                    PathManager.get_package_path())
        wopfile_path, wopfile_content = wopmars_runner.create_wopfile(path=wopfile_path)
        wopfile_content_bak = """rule SampleInformation:
    tool: vtam.wrapper.SampleInformation
    input:
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
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
        log_verbosity: 0


rule SortReads:
    tool: vtam.wrapper.SortReads
    input:
        table:
            Fasta: vtam.model.Fasta
            SampleInformation: vtam.model.SampleInformation
            Marker: vtam.model.Marker
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        file:
            sortreads: test/output/sortreads.tsv
    params:
        min_id: 0.8
        minseqlength: 32
        overhang: 0
        log_verbosity: 0


rule VariantReadCount:
    tool: vtam.wrapper.VariantReadCount
    input:
        file:
            sortreads: test/output/sortreads.tsv
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
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
        log_verbosity: 0


rule FilterLFN:
    tool: vtam.wrapper.FilterLFN
    input:
        table:
            Run: vtam.model.Run
            Marker: vtam.model.Marker
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            VariantReadCount: vtam.model.VariantReadCount
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
            threshold_specific: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            FilterLFN: vtam.model.FilterLFN
    params:
        filter_lfn_variant: 1
        lfn_variant_threshold: 0.001
        lfn_variant_replicate_threshold: 0.001
        lfn_biosample_replicate_threshold: 0.001
        lfn_read_count_threshold: 10
        log_verbosity: 0


rule FilterMinReplicateNumber:
    tool: vtam.wrapper.FilterMinReplicateNumber
    input:
        table:
            Run: vtam.model.Run
            Marker: vtam.model.Marker
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            FilterLFN: vtam.model.FilterLFN
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            FilterMinReplicateNumber: vtam.model.FilterMinReplicateNumber
    params:
        min_replicate_number: 2
        log_verbosity: 0


rule FilterPCRError:
    tool: vtam.wrapper.FilterPCRError
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Variant: vtam.model.Variant
            FilterMinReplicateNumber: vtam.model.FilterMinReplicateNumber
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            FilterPCRError: vtam.model.FilterPCRError
    params:
        pcr_error_var_prop: 0.1
        log_verbosity: 0


rule FilterChimera:
    tool: vtam.wrapper.FilterChimera
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Variant: vtam.model.Variant
            FilterPCRError: vtam.model.FilterPCRError
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            FilterChimera: vtam.model.FilterChimera
            FilterChimeraBorderline: vtam.model.FilterChimeraBorderline
    params:
        log_verbosity: 0


rule FilterRenkonen:
    tool: vtam.wrapper.FilterRenkonen
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            FilterChimera: vtam.model.FilterChimera
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            FilterRenkonen: vtam.model.FilterRenkonen
    params:
        renkonen_threshold: 0.1
        log_verbosity: 0


rule FilterIndel:
    tool: vtam.wrapper.FilterIndel
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Variant: vtam.model.Variant
            FilterRenkonen: vtam.model.FilterRenkonen
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            FilterIndel: vtam.model.FilterIndel
    params:
        log_verbosity: 0


rule FilterCodonStop:
    tool: vtam.wrapper.FilterCodonStop
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Variant: vtam.model.Variant
            FilterIndel: vtam.model.FilterIndel
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            FilterCodonStop: vtam.model.FilterCodonStop
    params:
        genetic_table_number: 5
        log_verbosity: 0


rule ReadCountAverageOverReplicates:
    tool: vtam.wrapper.ReadCountAverageOverReplicates
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            FilterCodonStop: vtam.model.FilterCodonStop
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            ReadCountAverageOverReplicates: vtam.model.ReadCountAverageOverReplicates
    params:
        log_verbosity: 0


rule TaxAssign:
    tool: vtam.wrapper.TaxAssign
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Variant: vtam.model.Variant
            FilterCodonStop: vtam.model.FilterCodonStop
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
            taxonomy: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        table:
            TaxAssign: vtam.model.TaxAssign
    params:
        identity_threshold: 97
        include_prop: 90
        min_number_of_taxa: 3
        blast_db: test/test_files/blastdb
        num_threads: 8
        log_verbosity: 0


rule MakeOtuTable:
    tool: vtam.wrapper.MakeOtuTable
    input:
        table:
            Marker: vtam.model.Marker
            Run: vtam.model.Run
            Biosample: vtam.model.Biosample
            Replicate: vtam.model.Replicate
            Variant: vtam.model.Variant
            FilterChimeraBorderline: vtam.model.FilterChimeraBorderline
            FilterCodonStop: vtam.model.FilterCodonStop
            TaxAssign: vtam.model.TaxAssign
        file:
            fastainfo: test/utils/test_wopmars_runner_wopfile_otu.py
            taxonomy: test/utils/test_wopmars_runner_wopfile_otu.py
    output:
        file:
            OTUTable: test/output/otutable.tsv
    params:
        log_verbosity: 0"""
        self.assertTrue(wopfile_content == wopfile_content_bak)
