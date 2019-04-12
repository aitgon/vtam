# -*- coding: utf-8 -*-
import os

import pandas
from unittest import TestCase

from numpy import nan

from wopmetabarcoding.utils.PathFinder import PathFinder
from wopmetabarcoding.wrapper.TaxAssignUtilities import f01_taxonomy_sqlite_to_df, f03_1_tax_id_to_taxonomy_lineage, \
    f05_select_ltg, \
    f03_import_blast_output_into_df, f04_blast_result_subset


class TestTaxAssign(TestCase):

    def setUp(self):
        self.taxonomy_db_df = f01_taxonomy_sqlite_to_df('tax.sqlite')
        #
        self.identity_threshold = 97
        #
        self.__testdir_path = os.path.join(PathFinder.get_module_test_path())
        self.blast_MFZR_002737_tsv = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "test_files", "blast_MFZR_002737.tsv")
        self.blast_MFZR_001274_tsv = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "test_files", "blast_MFZR_001274.tsv")
        self.v1_fasta = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "test_files", "MFZR_001274.fasta")

    def test_f03_import_blast_output_into_df(self):
        """This test assess whether the blast has been well imported into a df"""
        #
        blast_out_df = f03_import_blast_output_into_df(self.blast_MFZR_002737_tsv)
        #
        self.assertTrue(blast_out_df.to_dict('list') == {'target_id': [1049499563, 1049496963, 1049491687, 1049490545],
                                         'identity': [100.0, 100.0, 99.429, 99.429],
                                         'target_tax_id': [761875, 761875, 761875, 761875]})

    def test_f03_1_tax_id_to_taxonomy_lineage(self):
        tax_id = 183142
        #
        taxonomy_lineage_dic = f03_1_tax_id_to_taxonomy_lineage(tax_id, self.taxonomy_db_df)
        self.assertTrue({'tax_id': 183142, 'species': 183142, 'genus': 10194, 'family': 10193, 'order': 84394,
                         'superorder': 1709201, 'class': 10191, 'phylum': 10190, 'no rank': 131567, 'kingdom': 33208,
                         'superkingdom': 2759} == taxonomy_lineage_dic)

    def test_f04_blast_result_subset(self):
        # From
        # 'variant_id': 'MFZR_001274',
        # 'variant_sequence': 'TTTATACTTTATTTTTGGTGTTTGAGCCGGAATAATTGGCTTAAGAATAAGCCTGCTAATCCGTTTAGAGCTTGGGGTTCTATGACCCTTCCTAGGAGATGAGCATTTGTACAATGTCATCGTTACCGCTCATGCTTTTATCATAATTTTTTTTATGGTTATTCCAATTTCTATA',
        blast_out_dic = {
            'target_id': [514884684, 514884680],
            'identity': [80, 80],
            'target_tax_id': [1344033, 1344033]}
        blast_result_subset_df = pandas.DataFrame(data=blast_out_dic)
        tax_lineage_df = f04_blast_result_subset(blast_result_subset_df, self.taxonomy_db_df)
        self.assertTrue(tax_lineage_df.to_dict('list')=={'identity': [80, 80], 'class': [10191, 10191],
            'family': [204743, 204743], 'genus': [360692, 360692], 'kingdom': [33208, 33208],
            'no rank': [131567, 131567], 'order': [84394, 84394], 'phylum': [10190, 10190],
            'species': [1344033, 1344033], 'superkingdom': [2759, 2759], 'superorder': [1709201, 1709201], 'tax_id': [1344033, 1344033]})

    def test_f05_select_ltg_identity_80(self):
        # List of lineages that will correspond to list of tax_ids: One lineage per row
        tax_lineage_df = pandas.DataFrame(data={
            'species' : [666, 183142, 183142, 183142],
            'genus' : [10194, 10194, 10194, 10194],
            'order' : [10193, 10193, 10193, 10193],
            'superorder' : [84394, 84394, 84394, 84394],
        })
        identity = 80
        #
        ltg_tax_id, ltg_rank = f05_select_ltg(tax_lineage_df, identity)
        #
        self.assertTrue(ltg_tax_id == 183142)
        self.assertTrue(ltg_rank == 'species')

    def test_f05_select_ltg_identity_100(self):
        # List of lineages that will correspond to list of tax_ids: One lineage per row
        tax_lineage_df = pandas.DataFrame(data={
            'species' : [666, 183142, 183142, 183142],
            'genus' : [10194, 10194, 10194, 10194],
            'order' : [10193, 10193, 10193, 10193],
            'superorder' : [84394, 84394, 84394, 84394],
        })
        identity = 100
        #
        ltg_tax_id, ltg_rank = f05_select_ltg(tax_lineage_df, identity)
        #
        self.assertTrue(ltg_tax_id == 10194)
        self.assertTrue(ltg_rank == 'genus')

    def test_99_full_tax_assign_after_blast_MFZR_002737(self):
        """
        This test takes a blast result and return ltg_tax_id, ltg_rank at given identity
        """
        #
        blast_out_df = f03_import_blast_output_into_df(self.blast_MFZR_002737_tsv)
        #
        identity = 100
        blast_out_df.loc[blast_out_df.identity >= self.identity_threshold]
        blast_result_subset_df = blast_out_df.loc[blast_out_df.identity >= identity, ['target_id', 'target_tax_id']]
        tax_lineage_df = f04_blast_result_subset(blast_result_subset_df, self.taxonomy_db_df)
        ltg_tax_id, ltg_rank = f05_select_ltg(tax_lineage_df, identity=identity, identity_threshold=self.identity_threshold)
        self.assertTrue(ltg_tax_id==761875)
        self.assertTrue(ltg_rank=='species')

    def test_99_full_tax_assign_after_blast_MFZR_001274(self):
        """
        This test takes a blast result and return ltg_tax_id, ltg_rank at given identity
        """
        #
        blast_out_df = f03_import_blast_output_into_df(self.blast_MFZR_001274_tsv)
        #
        identity = 80
        blast_out_df.loc[blast_out_df.identity >= self.identity_threshold]
        blast_result_subset_df = blast_out_df.loc[blast_out_df.identity >= identity, ['target_id', 'target_tax_id']]
        tax_lineage_df = f04_blast_result_subset(blast_result_subset_df, self.taxonomy_db_df)
        ltg_tax_id, ltg_rank = f05_select_ltg(tax_lineage_df, identity=identity, identity_threshold=self.identity_threshold)
        self.assertTrue(ltg_tax_id==1344033)
        self.assertTrue(ltg_rank=='species')


    def test_f06_biopython_blast(self):
        from Bio.Blast import NCBIWWW
        with open(self.v1_fasta) as fin:
            sequence_data = fin.read()
            result_handle = NCBIWWW.qblast("blastn", "nt", sequence_data, ncbi_gi=True, format_type = 'Tabular')
            with open('results.tsv', 'w') as out_handle:
                # blast_v1_result = result_handle.read()
                out_handle.write(result_handle.read())
            result_handle.close()
        # https://github.com/darcyabjones/gi-to-tax/blob/master/README.md
        import pdb;pdb.set_trace()
