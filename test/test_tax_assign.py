# -*- coding: utf-8 -*-
import sqlite3
import pandas
from unittest import TestCase

from bin.create_taxonomy_db import create_parser, f_create_taxonomy_db

from wopmetabarcoding.utils.constants import rank_hierarchy


class TestCreateTaxonomyDBSqlite(TestCase):

    def setUp(self):
        self.parser = create_parser()

    def test_tax_id_to_taxonomy_lineage(self):
        # {'tax_id': 183142, 'species': 183142, 'genus': 10194, 'family': 10193}
        pass

    def test_tax_id_list_to_taxonomy_lineage_df(self):
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
        tax_id_list = list(set(blast_out_dic['tax_id']))

        tax_id_list1 = ( 1419335,1419335,1419335)

        lineage_dic = {

            'tax_id': [],

        }
        tax_id_list = list(set(blast_out_dic['tax_id']))

        # Todo: take 183142, and create dictionnary like:

        cnx = sqlite3.connect('tax.sqlite')
        df = pandas.read_sql_query("SELECT * FROM taxonomy", cnx)
        lineage_list = []
        tax_lineage_header = ['tax_seq_id'] + rank_hierarchy

        for tax_seq_id in tax_id_list:
            lineage_dic['tax_id'] = tax_seq_id
            while True:
                df_tax_id = df.loc[df.tax_id == tax_seq_id].copy()
                tax_id = df_tax_id['tax_id'].values[0]
                rank = df_tax_id['rank'].values[0]
                tax_seq_id = df_tax_id['parent_tax_id'].values[0]
                lineage_dic[rank] = tax_id
                # if tax_id equal 1 so that the last rank with his further tax_id
                if (tax_id == 1):
                    df_tax_id = df.loc[df.tax_id == tax_id].copy()
                    tax_id = df_tax_id['tax_id'].values[0]
                    rank = df_tax_id['rank'].values[0]
                    lineage_dic[rank] = tax_id
                    lineage_list.append(lineage_dic)
                    break
            # free al lvalue from our dictionnary
            lineage_dic = {
                'tax_id': [],
            }

        # convert the lineage_list to a lineage_df data frame

        lineage_df = pandas.DataFrame(lineage_list)

        import pdb;
        pdb.set_trace()

        # {'tax_id': 183142, 'species': 183142, 'genus': 10194, 'family': 10193, ''}
        #  Check code in TaxAssignUtilities:create_phylogenetic_line_df

        #
        self.assertTrue(lineage_df.loc[(lineage_df.tax_id == 183142),
                                       'genus'].values[0] == 10194)

        #
        self.assertTrue(lineage_df.loc[(lineage_df.tax_id == 183142),
                                       'family'].values[0] == 10193)
