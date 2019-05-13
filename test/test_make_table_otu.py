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


class TestMakeTableOTU(TestCase):

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
        # Initialize out_df
        otu_df = variant_df.copy()
        #
        # Add Variant Sequence length
        variant_df_tmp = variant_df.copy()
        variant_df_tmp['sequence_length'] = variant_df_tmp['variant_sequence'].str.len()
        otu_df = otu_df.merge(variant_df_tmp, on=['variant_id', 'variant_sequence'])
        #
        # Add read_count_sum_per_variant
        read_count_sum_per_variant = filter_codon_stop_df.groupby('variant_id').sum().reset_index()[
            ['variant_id', 'read_count']]
        otu_df = otu_df.merge(read_count_sum_per_variant, on='variant_id')
        #
        # Merge variants and read_count per biosample
        # otu_df = variant_df.merge(filter_codon_stop_df, on='variant_id')
        #
        ############################
        #
        # Prepare biosamples data
        #
        ############################
        #
        # Sum read counts over replicates of each biosample
        otu_biosamples_df = filter_codon_stop_df.groupby(['run_id', 'marker_id', 'variant_id', 'biosample_id'])[
            'read_count'].sum().reset_index()
        #
        # Replace biosample ids with their name
        otu_biosamples_df = otu_biosamples_df.merge(biosample_df, left_on='biosample_id', right_on='id')
        otu_biosamples_df.drop(['biosample_id', 'id'], axis=1, inplace=True)
        otu_biosamples_df = otu_biosamples_df.rename(columns={'name': 'biosample_name'})
        #
        # Replace marker ids with their name
        otu_biosamples_df = otu_biosamples_df.merge(marker_df, left_on='marker_id', right_on='id')
        otu_biosamples_df.drop(['marker_id', 'id'], axis=1, inplace=True)
        otu_biosamples_df = otu_biosamples_df.rename(columns={'name': 'marker_name'})
        #
        # Replace run ids with their name
        otu_biosamples_df = otu_biosamples_df.merge(run_df, left_on='run_id', right_on='id')
        otu_biosamples_df.drop(['run_id', 'id'], axis=1, inplace=True)
        otu_biosamples_df = otu_biosamples_df.rename(columns={'name': 'run_name'})
        #
        # Pivot biosamples
        otu_biosamples_df = otu_biosamples_df.pivot_table(index=['run_name', 'marker_name', 'variant_id'], columns="biosample_name",
                                        values='read_count').reset_index()
        #
        ############################
        #
        # Merge variant and biosample information
        #
        ############################
        otu_df = otu_df.merge(otu_biosamples_df, on='variant_id')
        #
        #
        # Merge ltg tax assign results
        # import pdb; pdb.set_trace()
        otu_df = otu_df.merge(ltg_tax_assign_df, on='variant_id')
        # getting the taxonomy_db to df
        con = sqlite3.connect(taxonomy_sqlite_path)
        sql = """SELECT *  FROM taxonomy """
        taxonomy_db_df = pandas.read_sql(sql=sql, con=con)
        con.close()

        ##########################
        #
        # Create lineage from each ltg_tax_id
        #
        ##########################
        list_lineage = []
        for tax_id in otu_df['ltg_tax_id'].tolist():
            dic_lineage = f04_1_tax_id_to_taxonomy_lineage(tax_id, taxonomy_db_df, give_tax_name=True)
            list_lineage.append(dic_lineage)
        lineage_df = pandas.DataFrame(data=list_lineage)
        lineage_list_df_columns_sorted = list(filter(lambda x: x in lineage_df.columns.tolist(), rank_hierarchy_otu_table))
        lineage_df = lineage_df[lineage_list_df_columns_sorted + ['tax_id']]


        # Merge Otu_df with the lineage df
        otu_df = otu_df.merge(lineage_df, left_on='ltg_tax_id', right_on='tax_id')
        otu_df.drop('tax_id', axis=1, inplace=True)

        # assert

        otu_df = f16_otu_table_maker(run_df, marker_df, variant_df, biosample_df, filter_codon_stop_df, ltg_tax_assign_df,taxonomy_db_df)


        self.assertTrue(otu_df.loc[(otu_df.variant_id == 15)
                                   & (otu_df.read_count == 120),
                                  'class'].values[0])=='Monogononta'

        self.assertTrue(otu_df.loc[(otu_df.variant_id == 15)
                                    & (otu_df.read_count == 120),
                                   'order'].values[0]) == 'Ploima'


