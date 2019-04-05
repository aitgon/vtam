# -*- coding: utf-8 -*-
import sqlite3
import pandas
from unittest import TestCase

from bin.create_taxonomy_db import create_parser, f_create_taxonomy_db

from wopmetabarcoding.utils.constants import rank_hierarchy
from wopmetabarcoding.wrapper.TaxAssignUtilities import taxonomy_sqlite_to_df, tax_id_to_taxonomy_lineage


class TestTaxAssign(TestCase):

    def setUp(self):
        self.taxonomy_db_df = taxonomy_sqlite_to_df('tax.sqlite')

    def test_tax_id_to_taxonomy_lineage(self):
        tax_id = 183142
        #
        taxonomy_lineage_dic = tax_id_to_taxonomy_lineage(tax_id, self.taxonomy_db_df)
        self.assertTrue({'tax_id': 183142, 'species': 183142, 'genus': 10194, 'family': 10193, 'order': 84394,
                         'superorder': 1709201, 'class': 10191, 'phylum': 10190, 'no rank': 131567, 'kingdom': 33208,
                         'superkingdom': 2759} == taxonomy_lineage_dic)

    def test_tax_id_list_to_taxonomy_lineage_df(self):
        #variant_id 1 from the emese file
        blast_out_dic = {
            'target_id': [187475693, 187475691, 187475699, 187475689, 347466873, 347466869, 347466867, 401880215,
                          401880085, 401880083, 1189131099, 1189131147, 1189131101, 401880223, 401880011, 514885218,
                          514885058, 514884984, 514884980, 514884850, 514884848, 514884842, 514884840, 514884834,
                          514884826, 514884824, 514884822, 514884820, 514884818, 514884774, 514884684, 514884680],
            'identity': [82.759, 82.759, 82.184, 82.184, 82.081, 82.081, 82.081, 81.714, 80.702, 80.702, 80.571, 80.46,
                         80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80],
            'tax_id': [183142, 183142, 183142, 183142, 1087463, 1087463, 1087463, 1213204, 1213227, 1213227, 1957050,
                       1957061, 1957050, 1225233, 1225250, 1344032, 1344033, 1344033, 1344033, 1344033, 1344033,
                       1344033, 1344033, 1344033, 1344033, 1344033, 1344033, 1344033, 1344033, 1344032, 1344033,
                       1344033]}
        #variant_id 3 from the emese file
        blast_out_dic2 = {
            'target_id': [1109411676, 1109411672, 1109411670,1109411663, 674660246 ],
            'identity': [100,100,100,100,100],
            'tax_id': [1419335,1419335,1419335,1419335,1419335]}

        blast_out_df = pandas.DataFrame(data=blast_out_dic2)
        #
        lineage_list = []
        for tax_id in set(blast_out_dic2['tax_id']):
            print(tax_id)
            lineage_list.append(tax_id_to_taxonomy_lineage(tax_id, self.taxonomy_db_df))
        lineage_df = pandas.DataFrame(data=lineage_list)
        lineage_list_df_columns = lineage_df.columns.tolist()
        rank_hierarchy_tax_id = ['tax_id'] + rank_hierarchy
        lineage_list_df_columns_sorted = list(filter(lambda x: x in lineage_list_df_columns, rank_hierarchy_tax_id))
        lineage_df = lineage_df[lineage_list_df_columns_sorted]
        variant_id_to_target_lineage_df = blast_out_df.merge(lineage_df, on='tax_id')



        #

        # variant_id_to_target_lineage_df['identity'].value_counts() / variant_id_to_target_lineage_df['identity'].count()
        import pdb;pdb.set_trace()

    # def test_tax_id_list_to_taxonomy_lineage_df_bak(self):
    #     blast_out_dic = {
    #         'target_id': [187475693, 187475691, 187475699, 187475689, 347466873, 347466869, 347466867, 401880215,
    #                       401880085, 401880083, 1189131099, 1189131147, 1189131101, 401880223, 401880011, 514885218,
    #                       514885058, 514884984, 514884980, 514884850, 514884848, 514884842, 514884840, 514884834,
    #                       514884826, 514884824, 514884822, 514884820, 514884818, 514884774, 514884684, 514884680],
    #         'identity': [82.759, 82.759, 82.184, 82.184, 82.081, 82.081, 82.081, 81.714, 80.702, 80.702, 80.571, 80.46,
    #                      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80],
    #         'tax_id': [183142, 183142, 183142, 183142, 1087463, 1087463, 1087463, 1213204, 1213227, 1213227, 1957050,
    #                    1957061, 1957050, 1225233, 1225250, 1344032, 1344033, 1344033, 1344033, 1344033, 1344033,
    #                    1344033, 1344033, 1344033, 1344033, 1344033, 1344033, 1344033, 1344033, 1344032, 1344033,
    #                    1344033]}
    #     tax_id_list = list(set(blast_out_dic['tax_id']))
    #
    #     tax_id_list1 = ( 1419335,1419335,1419335)
    #
    #     lineage_dic = {
    #
    #         'tax_id': [],
    #
    #     }
    #     tax_id_list = list(set(blast_out_dic['tax_id']))
    #
    #     # Todo: take 183142, and create dictionnary like:
    #
    #     cnx = sqlite3.connect('tax.sqlite')
    #     df = pandas.read_sql_query("SELECT * FROM taxonomy", cnx)
    #     lineage_list = []
    #     tax_lineage_header = ['tax_seq_id'] + rank_hierarchy
    #     for tax_seq_id in tax_id_list:
    #         lineage_dic['tax_id'] = tax_seq_id
    #         while True:
    #             df_tax_id = df.loc[df.tax_id == tax_seq_id].copy()
    #             tax_id = df_tax_id['tax_id'].values[0]
    #             rank = df_tax_id['rank'].values[0]
    #             tax_seq_id = df_tax_id['parent_tax_id'].values[0]
    #             lineage_dic[rank] = tax_id
    #             # if tax_id equal 1 so that the last rank with his further tax_id
    #             if (tax_id == 1):
    #                 df_tax_id = df.loc[df.tax_id == tax_id].copy()
    #                 tax_id = df_tax_id['tax_id'].values[0]
    #                 rank = df_tax_id['rank'].values[0]
    #                 lineage_dic[rank] = tax_id
    #                 lineage_list.append(lineage_dic)
    #                 break
    #         # free al lvalue from our dictionnary
    #         lineage_dic = {
    #             'tax_id': [],
    #         }
    #         import pdb; pdb.set_trace()
    #     # convert the lineage_list to a lineage_df data frame
    #
    #     lineage_df = pandas.DataFrame(lineage_list)
    #
    #     import pdb;
    #     pdb.set_trace()
    #
    #     # {'tax_id': 183142, 'species': 183142, 'genus': 10194, 'family': 10193, ''}
    #     #  Check code in TaxAssignUtilities:create_phylogenetic_line_df
    #
    #     #
    #     self.assertTrue(lineage_df.loc[(lineage_df.tax_id == 183142),
    #                                    'genus'].values[0] == 10194)
    #
    #     #
    #     self.assertTrue(lineage_df.loc[(lineage_df.tax_id == 183142),
    #                                    'family'].values[0] == 10193)
