import inspect
import numpy
import os
import pandas
import pathlib
import unittest

from vtam.utils.PathManager import PathManager
from vtam.CommandTaxonomy import CommandTaxonomy
from vtam.utils.Logger import Logger
from vtam.utils.VariantDF import VariantDF
from vtam.utils.TaxAssignRunner import f06_select_ltg


class TestTaxAssign(unittest.TestCase):

    def setUp(self):
        #####################################
        #
        # Download taxonomy.tsv
        #
        #####################################
        #
        self.ltg_rule_threshold = 97
        self.min_number_of_taxa = 3
        self.include_prop = 90
        #
        self.__testdir_path = PathManager.get_test_path()

    @classmethod
    def setUpClass(cls):
        """ get_some_resource() is slow, to avoid calling it for each tests use setUpClass()
            and store the result as class variable
        """
        super(TestTaxAssign, cls).setUpClass()
        # create_vtam_data_dir()
        testdir_path = os.path.join(PathManager.get_test_path())
        cls.outdir_path = os.path.join(testdir_path, "outdir")
        pathlib.Path(cls.outdir_path).mkdir(exist_ok=True, parents=True)
        taxonomy_tsv_path = os.path.join(cls.outdir_path, "taxonomy.tsv")
        CommandTaxonomy(
            taxonomy_tsv=taxonomy_tsv_path).download_precomputed_taxonomy()

    def test_variant_df_to_fasta(self):
        variant_dic = {
            'id': [
                57,
                107],
            'sequence': [
                'ACTATATTTTATTTTTGGGGCTTGATCCGGAATGCTGGGCACCTCTCTAAGCCTTCTAATTCGTGCCGAGCTGGGGCACCCGGGTTCTTTAATTGGCGACGATCAAATTTACAATGTAATCGTCACAGCCCATGCTTTTATTATGATTTTTTTCATGGTTATGCCTATTATAATC',
                'ACTTTATTTCATTTTCGGAACATTTGCAGGAGTTGTAGGAACTTTACTTTCATTATTTATTCGTCTTGAATTAGCTTATCCAGGAAATCAATTTTTTTTAGGAAATCACCAACTTTATAATGTGGTTGTGACAGCACATGCTTTTATCATGATTTTTTTCATGGTTATGCCGATTTTAATC']}
        variant_df = pandas.DataFrame(data=variant_dic)
        variant_df.set_index('id', inplace=True)

        #
        Logger.instance().debug(
            "file: {}; line: {}; Create SortedReadFile from Variants".format(
                __file__, inspect.currentframe().f_lineno, 'PoolMarkers'))
        this_tempdir = os.path.join(
            PathManager.instance().get_tempdir(),
            os.path.basename(__file__))
        pathlib.Path(this_tempdir).mkdir(exist_ok=True)
        variant_fasta = os.path.join(this_tempdir, 'variant.fasta')
        variant_df_utils = VariantDF(variant_df)
        variant_df_utils.to_fasta(fasta_path=variant_fasta)
        #
        variant_fasta_content_expected = """>57\nACTATATTTTATTTTTGGGGCTTGATCCGGAATGCTGGGCACCTCTCTAAGCCTTCTAATTCGTGCCGAGCTGGGGCACCCGGGTTCTTTAATTGGCGACGATCAAATTTACAATGTAATCGTCACAGCCCATGCTTTTATTATGATTTTTTTCATGGTTATGCCTATTATAATC\n>107\nACTTTATTTCATTTTCGGAACATTTGCAGGAGTTGTAGGAACTTTACTTTCATTATTTATTCGTCTTGAATTAGCTTATCCAGGAAATCAATTTTTTTTAGGAAATCACCAACTTTATAATGTGGTTGTGACAGCACATGCTTTTATCATGATTTTTTTCATGGTTATGCCGATTTTAATC\n"""
        with open(variant_fasta, 'r') as fin:
            variant_fasta_content = fin.read()
        self.assertTrue(variant_fasta_content_expected ==
                        variant_fasta_content)

    def test_f06_select_ltg_identity_80(self):
        # List of lineages that will correspond to list of tax_ids: One lineage
        # per row
        tax_lineage_df = pandas.DataFrame(data={
            'species': [666, 183142, 183142, 183142],
            'genus': [10194, 10194, 10194, 10194],
            'order': [10193, 10193, 10193, 10193],
            'superorder': [84394, 84394, 84394, 84394],
        })
        identity = 80
        ltg_tax_id, ltg_rank = f06_select_ltg(
            tax_lineage_df=tax_lineage_df, include_prop=self.include_prop)
        #
        # import pdb; pdb.set_trace()
        self.assertTrue(ltg_tax_id == 10194)
        self.assertTrue(ltg_rank == 'genus')

    def test_f06_select_ltg_column_none(self):
        # List of lineages that will correspond to list of tax_ids: One lineage
        # per row
        tax_lineage_df = pandas.DataFrame(data={
            'species': [666, 183142, 183142, 183142],
            'subgenus': [numpy.nan] * 4,
            'genus': [10194, 10194, 10194, 10194],
            'order': [10193, 10193, 10193, 10193],
            'superorder': [84394, 84394, 84394, 84394],
        })
        #
        identity = 80
        #
        ltg_tax_id, ltg_rank = f06_select_ltg(
            tax_lineage_df=tax_lineage_df, include_prop=self.include_prop)
        #
        self.assertTrue(ltg_tax_id == 10194)
        self.assertTrue(ltg_rank == 'genus')

    def test_f05_select_ltg_identity_100(self):
        # List of lineages that will correspond to list of tax_ids: One lineage
        # per row
        tax_lineage_df = pandas.DataFrame(data={
            'species': [666, 183142, 183142, 183142],
            'genus': [10194, 10194, 10194, 10194],
            'order': [10193, 10193, 10193, 10193],
            'superorder': [84394, 84394, 84394, 84394],
        })
        identity = 100
        #
        ltg_tax_id, ltg_rank = f06_select_ltg(
            tax_lineage_df=tax_lineage_df, include_prop=self.include_prop)
        #
        self.assertTrue(ltg_tax_id == 10194)
        self.assertTrue(ltg_rank == 'genus')
