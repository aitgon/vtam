# -*- coding: utf-8 -*-
import inspect
import os
import sqlite3

from unittest import TestCase

import pandas

from wopmetabarcoding.utils.constants import rank_hierarchy_otu_table
from wopmetabarcoding.wrapper.TaxAssignUtilities import f04_1_tax_id_to_taxonomy_lineage

from wopmetabarcoding.wrapper.TaxAssignUtilities import f01_taxonomy_sqlite_to_df

from wopmetabarcoding.wrapper.MakeOtuTable import f16_otu_table_maker


class TestMakeOtuTable(TestCase):

    def test_f01_make_table_out(self):
        #  Input
        try:
            taxonomy_sqlite_path = os.path.join(os.environ['DIR_DATA_NON_GIT'], 'taxonomy.sqlite')
        except:
            raise Exception('Please, set the DIR_DATA_NON_GIT environment variable. See the manual')

        run_dic = {
            'id': [1],
            'name': ['prerun']
        }
        marker_dic = {
            'id': [1],
            'name': ['MFZR']
        }
        variant_dic = {
            'variant_id': [30,15],
            'variant_sequence': [
                'ACTATATTTTATTTTTGGGGCTTGATCCGGAATGCTGGGCACCTCTCTAAGCCTTCTAATTCGTGCCGAGCTGGGGCACCCGGGTTCTTTAATTGGCGACGATCAAATTTACAATGTAATCGTCACAGCCCATGCTTTTATTATGATTTTTTTCATGGTTATGCCTATTATAATC',
                'TTTATACTTTATTTTTGGTGTTTGAGCCGGAATAATTGGCTTAAGAATAAGCCTGCTAATCCGTTTAGAGCTTGGGGTTCTATGACCCTTCCTAGGAGATGAGCATTTGTACAATGTCATCGTTACCGCTCATGCTTTTATCATAATTTTTTTTATGGTTATTCCAATTTCTATA']
        }
        filter_codon_stop_dic = {'run_id': [1, 1, 1, 1, 1, 1,1,1,1], 'marker_id': [1, 1, 1, 1, 1, 1,1,1,1],
                                 'variant_id': [30, 30, 30, 30, 30, 30,15,15,15], 'biosample_id': [1, 1, 1, 2, 2, 2,1,2,1],
                                 'replicate_id': [1, 2, 3, 1, 2, 3,1,2,3], 'read_count': [183, 93, 42, 175, 31, 63,20,40,60],
                                 'fiter_id': [14, 14, 14, 14, 14, 14,14,14,14], 'filter_delete': [0, 0, 0, 0, 0, 0,0,0,0]}
        #
        biosample_dic = {'id': [1, 2], 'name': ['Tpos1_prerun', 'Tpos2_prerun']}
        #
        ltg_tax_assign_dic = {'variant_id': [30,15], 'identity': [85,80], 'ltg_rank': ['species','species'], 'ltg_tax_id': [268290,84394]}

        #
        #  Get tables/df
        run_df = pandas.DataFrame(run_dic, index=None)
        marker_df = pandas.DataFrame(marker_dic, index=None)
        variant_df = pandas.DataFrame(variant_dic, index=None)
        biosample_df = pandas.DataFrame(biosample_dic, index=None)
        filter_codon_stop_df = pandas.DataFrame(filter_codon_stop_dic, index=None)
        ltg_tax_assign_df = pandas.DataFrame(ltg_tax_assign_dic, index=None)
        #
        #
        # getting the taxonomy_db to df
        con = sqlite3.connect(taxonomy_sqlite_path)
        sql = """SELECT *  FROM taxonomy """
        taxonomy_db_df = pandas.read_sql(sql=sql, con=con)
        con.close()

        otu_df = f16_otu_table_maker(run_df, marker_df, variant_df, biosample_df, filter_codon_stop_df, ltg_tax_assign_df, taxonomy_db_df)

        # otu_final_df = otu_df[columns]
        self.assertTrue(otu_df.loc[(otu_df.variant_id == 15)
                                   & (otu_df.read_count == 120),
                                  'class'].values[0])=='Monogononta'

        self.assertTrue(otu_df.loc[(otu_df.variant_id == 15)
                                    & (otu_df.read_count == 120),
                                   'order'].values[0]) == 'Ploima'
