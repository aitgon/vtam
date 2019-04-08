# -*- coding: utf-8 -*-

import pandas
from unittest import TestCase
from wopmetabarcoding.wrapper.TaxAssignUtilities import taxonomy_sqlite_to_df, tax_id_to_taxonomy_lineage, select_ltg


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

    def test_select_ltg_identity_80(self):
        # List of lineages that will correspond to list of tax_ids: One lineage per row
        tax_lineage_df = pandas.DataFrame(data={
            'species' : [666, 183142, 183142, 183142],
            'genus' : [10194, 10194, 10194, 10194],
            'order' : [10193, 10193, 10193, 10193],
            'superorder' : [84394, 84394, 84394, 84394],
        })
        identity = 80
        #
        ltg_tax_id, ltg_rank = select_ltg(tax_lineage_df, identity)
        #
        self.assertTrue(ltg_tax_id == 183142)
        self.assertTrue(ltg_rank == 'species')

    def test_select_ltg_identity_100(self):
        # List of lineages that will correspond to list of tax_ids: One lineage per row
        tax_lineage_df = pandas.DataFrame(data={
            'species' : [666, 183142, 183142, 183142],
            'genus' : [10194, 10194, 10194, 10194],
            'order' : [10193, 10193, 10193, 10193],
            'superorder' : [84394, 84394, 84394, 84394],
        })
        identity = 100
        #
        ltg_tax_id, ltg_rank = select_ltg(tax_lineage_df, identity)
        #
        self.assertTrue(ltg_tax_id == 10194)
        self.assertTrue(ltg_rank == 'genus')

    def test_tax_id_list_to_taxonomy_lineage_df(self):
        #
        # 'variant_id': 'MFZR_001274',
        # 'variant_sequence': 'TTTATACTTTATTTTTGGTGTTTGAGCCGGAATAATTGGCTTAAGAATAAGCCTGCTAATCCGTTTAGAGCTTGGGGTTCTATGACCCTTCCTAGGAGATGAGCATTTGTACAATGTCATCGTTACCGCTCATGCTTTTATCATAATTTTTTTTATGGTTATTCCAATTTCTATA',
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
        blast_out_df = pandas.DataFrame(data=blast_out_dic)
        # l = pandas.Series(blast_out_df.tax_id.unique()).apply(lambda x: tax_id_to_taxonomy_lineage(x, self.taxonomy_db_df))
        lineage_list = []
        for tax_id in blast_out_df.tax_id.unique().tolist():
            lineage_list.append(tax_id_to_taxonomy_lineage(tax_id, self.taxonomy_db_df))
        tax_lineage_df=pandas.DataFrame(lineage_list)
        tax_lineage_df.index=tax_lineage_df.tax_id
        tax_lineage_df.drop('tax_id', axis=1, inplace=True)
        ltg_tax_id, ltg_rank = select_ltg(tax_lineage_df, identity=80)
        #
        self.assertTrue(ltg_tax_id == 204743,)
        self.assertTrue(ltg_rank == 'family')

    def test_tax_id_list_to_taxonomy_lineage_df_bak(self):
        #
        # 'variant_id': 'MFZR_001274',
        # 'variant_sequence': 'TTTATACTTTATTTTTGGTGTTTGAGCCGGAATAATTGGCTTAAGAATAAGCCTGCTAATCCGTTTAGAGCTTGGGGTTCTATGACCCTTCCTAGGAGATGAGCATTTGTACAATGTCATCGTTACCGCTCATGCTTTTATCATAATTTTTTTTATGGTTATTCCAATTTCTATA',
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
        blast_out_dic2 = {
            'variant_id': 'MFZR_001790',
            'variant_sequence': 'CTTATATTTTATTTTCGGAGCTTGAGCTGGAATAGTTGGGACTTCCCTTAGTATACTAATTCGAGCTGAATTAGGTCATCCTGGTTCACTTATTGGTGATGATCAAATTTATAATGTAATTGTTACAGCTCATGCTTTTGTAATAATTTTTTTTATAGTTATACCTATTATAATT',
            'target_id': [1109411676, 1109411672, 1109411670,1109411663, 674660246 ],
            'identity': [100,100,100,100,100],
            'tax_id': [1419335,1419335,1419335,1419335,1419335]}

        # variant_id 5 from the emese file
        blast_out_dic3 = {
            'variant_id': 'MFZR_002045',
            'variant_sequence': 'TTTATACTTTGTGTTTGGAGCTTGGGCTGGAATAGTCGGCTCCTCTCTAAGGGTTTTAATTCGTCTAGAGTTAGGGCAACCTGGCTCATTAATTGGGGATGATCAGATCTACAATGTAGTTGTGACAGCTCATGCTTTCGTTATAATCTTCTTTATGGTGATACCTGCTATAATT',
            'target_id': [1258178191, 1214941590, 1109411670 , 1214941112, 1214940769, 1214940743, 1214940345, 1214940329],
            'identity': [99.429, 99.429, 99.429, 99.429, 99.429, 99.429, 99.429, 99.429 ],
            'tax_id': [189839, 189839, 189839, 189839, 189839, 189839, 189839, 189839]}

        #var id 8 from emese file
        blast_out_dic4 = {
            'variant_id': 'MFZR_002737',
            'variant_sequence': 'TCTTTATTTTATATTTGGAATCTGATCAGGCTTAGTAGGATCATCTTTAAGTTTTATTATTCGAATAGAATTAAGAACTCCTGGAAGATTTATTGGAAATGACCAAATTTATAATGTTGTAGTAACTTCTCATGCTTTTATCATAATTTTTTTTATAGTAATGCCTATTATAATC',
            'target_id': [1049499563, 1049496963, 1239402664, 1239402658, 1049505397],
            'identity': [100, 100, 100, 100,100],
            'tax_id': [761875, 761875, 761875, 761875, 761875]}


