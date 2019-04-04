# -*- coding: utf-8 -*-

from unittest import TestCase

from bin.create_taxonomy_db import create_parser, f_create_taxonomy_db

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
        import pdb; pdb.set_trace()
        # Todo: take 183142, and create dictionnary like:
        # {'tax_id': 183142, 'species': 183142, 'genus': 10194, 'family': 10193, ''}
        # Check code in TaxAssignUtilities:create_phylogenetic_line_df
