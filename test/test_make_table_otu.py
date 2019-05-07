# -*- coding: utf-8 -*-
import inspect
import os
import sqlite3

from unittest import TestCase

import pandas

from wopmetabarcoding.wrapper.TaxAssignUtilities import f04_1_tax_id_to_taxonomy_lineage

from wopmetabarcoding.wrapper.TaxAssignUtilities import f01_taxonomy_sqlite_to_df


class TestMakeTableOTU(TestCase):

    def test_f01_make_table_out(self):

        # Input

        taxonomy_sqlite_path="taxonomy.sqlite"

        variant_dic = {
            'variant_id' : [30],
            'variant_sequence' : ['ACTATATTTTATTTTTGGGGCTTGATCCGGAATGCTGGGCACCTCTCTAAGCCTTCTAATTCGTGCCGAGCTGGGGCACCCGGGTTCTTTAATTGGCGACGATCAAATTTACAATGTAATCGTCACAGCCCATGCTTTTATTATGATTTTTTTCATGGTTATGCCTATTATAATC']
        }
        filter_codon_stop_dic = {'run_id': [1, 1, 1, 1, 1, 1], 'marker_id': [1, 1, 1, 1, 1, 1],
                                 'variant_id': [30, 30, 30, 30, 30, 30], 'biosample_id': [1, 1, 1, 2, 2, 2],
                                 'replicate_id': [1, 2, 3, 1, 2, 3], 'read_count': [183, 93, 42, 175, 31, 63],
                                 'fiter_id': [14, 14, 14, 14, 14, 14], 'filter_delete': [0, 0, 0, 0, 0, 0]}
        #
        biosample_dic = { 'id': [1, 2],'name': ['Tpos1_prerun','Tpos2_prerun']}
        #
        taxonomy_dic = { 'variant_id' : [30], 'identity': [85], 'ltg_rank': ['species'],  'ltg_tax_id':[268290]}

        #
        # Get tables
        variant_df = pandas.DataFrame(variant_dic, index=None)
        biosample_df = pandas.DataFrame(biosample_dic, index=None)
        filter_codon_stop_df = pandas.DataFrame(filter_codon_stop_dic, index=None)
        taxonomy_df = pandas.DataFrame(taxonomy_dic, index=None)
        # out = pandas.DataFrame(columns=['variant_df', 'sequence_length'])
        #
        # Merge table
        otu_df = variant_df.merge(filter_codon_stop_df, on='variant_id')
        #
        # Sum read counts over replicates of each biosample
        otu_df = otu_df.groupby(['run_id', 'marker_id', 'variant_id', 'variant_sequence', 'biosample_id'])[
            'read_count'].sum().reset_index()
        #
        # Pivot biosamples
        otu_tmp_df = otu_df.pivot_table(index=['run_id', 'marker_id', 'variant_id'], columns="biosample_id", values='read_count').reset_index()
        #
        # Merge pivoted biosamples with variant information
        otu_df = otu_df[['run_id', 'marker_id', 'variant_id', 'variant_sequence']].drop_duplicates().merge(otu_tmp_df,
                                                                            on=['run_id', 'marker_id', 'variant_id'])
        import pdb; pdb.set_trace()
        # # variant_id_delete_list = []
        # out_dic = {}
        # for variant_id in  variant_df["variant_id"].tolist():
        #     #sequence length
        #     out_dic['variant_id'] = variant_id
        #     df = variant_df.loc[variant_df.variant_id == variant_id]
        #     out_dic['sequence_length'] = df["variant_sequence"].str.len()
        #
        #     # out_dic['sample']
        #
        #     # out_dic[biosample_df]= biosample_df["name"]
        #
        #
        #     # ltg lineage - Ltg taxid - Ltg identity - LTG rank
        #     ltg_df = taxonomy_df.loc[taxonomy_df.variant_id == variant_id]
        #
        #     out_dic['ltg_tax_id'] = ltg_df["ltg_tax_id"]
        #     out_dic['ltg_rank'] = ltg_df["ltg_rank"]
        #     out_dic['identity'] = ltg_df["identity"]
        #
        #
        #     # taxonomy db
        #     # getting the taxonomy_db to df
        #     con = sqlite3.connect(taxonomy_sqlite_path)
        #     sql = """SELECT *  FROM taxonomy """
        #     taxonomy_db_df = pandas.read_sql(sql=sql, con=con)
        #     con.close()
        #     ltg_taxid = ltg_df["ltg_tax_id"]
        #
        #
        #     # dic_lineage = f04_1_tax_id_to_taxonomy_lineage(ltg_taxid, taxonomy_db_df)
        #     import pdb;
        #     pdb.set_trace()
        #
        #     # add variant-sequence
        #     df_sequence = variant_df.loc[variant_df.variant_id == variant_id]
        #     out_dic['variant_sequence'] = df_sequence["variant_sequence"]
        #
        #
        # for biosample_id in  biosample_df["id"].tolist():
        #     df_biosample =  biosample_df.loc[biosample_df.id ==biosample_id]
        #     # sample and read count
        #     vatriant_filter_codon_stop_df = filter_codon_stop_df.loc[filter_codon_stop_df.variant_id == variant_id]
        #     read_count_per_biosample = vatriant_filter_codon_stop_df.groupby(['biosample_id'])[
        #         'read_count'].sum().reset_index()
        #     df_biosample_read_count =  read_count_per_biosample.loc[read_count_per_biosample["biosample_id"] == biosample_id]
        #
        #     out_dic[df_biosample["name"]]= df_biosample_read_count["read_count"]
        #
        #
        # out_df = pandas.DataFrame(out_dic, index=None)
        # import pdb;
        # pdb.set_trace()