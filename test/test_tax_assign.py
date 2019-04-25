# -*- coding: utf-8 -*-
import errno
import inspect
import os
import sqlite3
import tarfile
import urllib

import numpy
import pandas
from unittest import TestCase

from Bio.Blast import NCBIWWW

from wopmetabarcoding.utils.PathFinder import PathFinder
from wopmetabarcoding.utils.constants import tempdir, data_dir, public_data_dir
from wopmetabarcoding.utils.logger import logger
from wopmetabarcoding.wrapper.TaxAssignUtilities import f01_taxonomy_sqlite_to_df, f04_1_tax_id_to_taxonomy_lineage, \
    f06_select_ltg, \
    f04_import_blast_output_into_df, f05_blast_result_subset, f02_variant_df_to_fasta

import gzip
import shutil


class TestTaxAssign(TestCase):

    def setUp(self):
        #####################################
        #
        # Download taxonomy.sqlite
        #
        #####################################

        file_remote = os.path.join(public_data_dir, "taxonomy.sqlite")
        taxonomy_sqlite_path = os.path.join(data_dir, os.path.basename(file_remote))
        if not os.path.isfile(taxonomy_sqlite_path):
            logger.debug(
                "file: {}; line: {}; Downloading taxonomy.sqlite".format(__file__, inspect.currentframe().f_lineno))
            urllib.request.urlretrieve(file_remote, taxonomy_sqlite_path)
        #
        self.taxonomy_db_df = f01_taxonomy_sqlite_to_df(taxonomy_sqlite_path)
        #
        self.identity_threshold = 97
        #
        self.__testdir_path = os.path.join(PathFinder.get_module_test_path())
        self.blast_MFZR_002737_tsv = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "test_files", "blast_MFZR_002737.tsv")
        self.blast_MFZR_001274_tsv = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "test_files", "blast_MFZR_001274.tsv")
        self.v1_fasta = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "MFZR_001274.fasta")
        self.blast_MFZR_002737_tsv = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "test_files",
                                                  "blast_MFZR_002737.tsv")

    def test_f02_variant_df_to_fasta(self):
        variant_dic = {
            'variant_id' : [57, 107],
            'variant_sequence' : ['ACTATATTTTATTTTTGGGGCTTGATCCGGAATGCTGGGCACCTCTCTAAGCCTTCTAATTCGTGCCGAGCTGGGGCACCCGGGTTCTTTAATTGGCGACGATCAAATTTACAATGTAATCGTCACAGCCCATGCTTTTATTATGATTTTTTTCATGGTTATGCCTATTATAATC'
                                  , 'ACTTTATTTCATTTTCGGAACATTTGCAGGAGTTGTAGGAACTTTACTTTCATTATTTATTCGTCTTGAATTAGCTTATCCAGGAAATCAATTTTTTTTAGGAAATCACCAACTTTATAATGTGGTTGTGACAGCACATGCTTTTATCATGATTTTTTTCATGGTTATGCCGATTTTAATC']
        }
        variant_df = pandas.DataFrame(data=variant_dic)
        #
        logger.debug(
            "file: {}; line: {}; Create Fasta from Variants".format(__file__, inspect.currentframe().f_lineno ,'TaxAssign'))
        this_tempdir = os.path.join(tempdir, os.path.basename(__file__))
        try:
            os.makedirs(this_tempdir)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        variant_fasta = os.path.join(this_tempdir, 'variant.fasta')
        f02_variant_df_to_fasta(variant_df, variant_fasta)
        #
        variant_fasta_content_expected =""">57\nACTATATTTTATTTTTGGGGCTTGATCCGGAATGCTGGGCACCTCTCTAAGCCTTCTAATTCGTGCCGAGCTGGGGCACCCGGGTTCTTTAATTGGCGACGATCAAATTTACAATGTAATCGTCACAGCCCATGCTTTTATTATGATTTTTTTCATGGTTATGCCTATTATAATC\n>107\nACTTTATTTCATTTTCGGAACATTTGCAGGAGTTGTAGGAACTTTACTTTCATTATTTATTCGTCTTGAATTAGCTTATCCAGGAAATCAATTTTTTTTAGGAAATCACCAACTTTATAATGTGGTTGTGACAGCACATGCTTTTATCATGATTTTTTTCATGGTTATGCCGATTTTAATC\n"""
        with open(variant_fasta, 'r') as fin:
            variant_fasta_content = fin.read()
        self.assertTrue(variant_fasta_content_expected == variant_fasta_content)

    #
    # Commented because too slow
    # Uncomment to test qblast
    # def test_f06_2_run_qblast(self):
    #     #
    #     # Run and read blast result
    #     with open(self.v1_fasta) as fin:
    #         variant_fasta_content = fin.read()
    #         logger.debug(
    #             "file: {}; line: {}; Blasting...".format(__file__, inspect.currentframe().f_lineno))
    #         result_handle = NCBIWWW.qblast("blastn", "nt", variant_fasta_content, format_type = 'Tabular')
    #         blast_result_tsv = os.path.join(tempdir, "tax_assign_blast.tsv")
    #         with open(blast_result_tsv, 'w') as out_handle:
    #             out_handle.write(result_handle.read())
    #         result_handle.close()
    #     blast_result_df = pandas.read_csv(blast_result_tsv, sep="\t", skiprows=13, usecols=[0, 1, 2],
    #                          header=None, names=['variant_id', 'gb_accession', 'identity'])

    def test_f06_3_annotate_blast_output_with_tax_id(self):
        #
        qblast_MFZR_001274_tsv = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "test_files",
                                                  "qblast_MFZR_001274.tsv")
        nucl_gb_accession2taxid_MFZR_00001274_sqlite = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "test_files",
                                                  "nucl_gb_accession2taxid_MFZR_001274.sqlite")
        #
        # Run and read blast result
        blast_result_df = pandas.read_csv(qblast_MFZR_001274_tsv, sep="\t", skiprows=13, usecols=[0, 1, 2],
                             header=None, names=['variant_id', 'gb_accession', 'identity'])
        blast_result_df = blast_result_df.loc[~pandas.isnull(blast_result_df).any(axis=1)]
        # import pdb; pdb.set_trace()
        con = sqlite3.connect(nucl_gb_accession2taxid_MFZR_00001274_sqlite)
        sql = """SELECT gb_accession, tax_id FROM nucl_gb_accession2taxid WHERE gb_accession IN {}""".format(tuple(blast_result_df.gb_accession.tolist()))
        gb_accession_to_tax_id_df = pandas.read_sql(sql=sql, con=con)
        con.close()
        #
        blast_result_tax_id_df = blast_result_df.merge(gb_accession_to_tax_id_df, on='gb_accession')


    def test_f03_import_blast_output_into_df(self):
        """This test assess whether the blast has been well imported into a df"""
        #
        blast_out_df = f04_import_blast_output_into_df(self.blast_MFZR_002737_tsv)
        #
        self.assertTrue(blast_out_df.to_dict('list') == {'target_id': [1049499563, 1049496963, 1049491687, 1049490545],
                                         'identity': [100.0, 100.0, 99.429, 99.429],
                                         'target_tax_id': [761875, 761875, 761875, 761875]})

    def test_f03_1_tax_id_to_taxonomy_lineage(self):
        tax_id = 183142
        #
        taxonomy_lineage_dic = f04_1_tax_id_to_taxonomy_lineage(tax_id, self.taxonomy_db_df)
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
        tax_lineage_df = f05_blast_result_subset(blast_result_subset_df, self.taxonomy_db_df)
        self.assertTrue(tax_lineage_df.to_dict('list')=={'identity': [80, 80], 'class': [10191, 10191],
            'family': [204743, 204743], 'genus': [360692, 360692], 'kingdom': [33208, 33208],
            'no rank': [131567, 131567], 'order': [84394, 84394], 'phylum': [10190, 10190],
            'species': [1344033, 1344033], 'superkingdom': [2759, 2759], 'superorder': [1709201, 1709201], 'tax_id': [1344033, 1344033]})

    def test_f06_select_ltg_identity_80(self):
        # List of lineages that will correspond to list of tax_ids: One lineage per row
        tax_lineage_df = pandas.DataFrame(data={
            'species' : [666, 183142, 183142, 183142],
            'genus' : [10194, 10194, 10194, 10194],
            'order' : [10193, 10193, 10193, 10193],
            'superorder' : [84394, 84394, 84394, 84394],
        })
        identity = 80
        #
        ltg_tax_id, ltg_rank = f06_select_ltg(tax_lineage_df, identity)
        #
        self.assertTrue(ltg_tax_id == 183142)
        self.assertTrue(ltg_rank == 'species')

    def test_f06_select_ltg_column_none(self):
        # List of lineages that will correspond to list of tax_ids: One lineage per row
        tax_lineage_df = pandas.DataFrame(data={
            'species' : [666, 183142, 183142, 183142],
            'subgenus' : [numpy.nan] * 4,
            'genus' : [10194, 10194, 10194, 10194],
            'order' : [10193, 10193, 10193, 10193],
            'superorder' : [84394, 84394, 84394, 84394],
        })
        #
        identity = 80
        #
        ltg_tax_id, ltg_rank = f06_select_ltg(tax_lineage_df, identity)
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
        ltg_tax_id, ltg_rank = f06_select_ltg(tax_lineage_df, identity)
        #
        self.assertTrue(ltg_tax_id == 10194)
        self.assertTrue(ltg_rank == 'genus')

    def test_99_full_tax_assign_after_blast_MFZR_002737(self):
        """
        This test takes a blast result and return ltg_tax_id, ltg_rank at given identity
        """
        #
        blast_out_df = f04_import_blast_output_into_df(self.blast_MFZR_002737_tsv)
        #
        identity = 100
        blast_out_df.loc[blast_out_df.identity >= self.identity_threshold]
        blast_result_subset_df = blast_out_df.loc[blast_out_df.identity >= identity, ['target_id', 'target_tax_id']]
        tax_lineage_df = f05_blast_result_subset(blast_result_subset_df, self.taxonomy_db_df)
        ltg_tax_id, ltg_rank = f06_select_ltg(tax_lineage_df, identity=identity, identity_threshold=self.identity_threshold)
        self.assertTrue(ltg_tax_id==761875)
        self.assertTrue(ltg_rank=='species')

    def test_99_full_tax_assign_after_blast_MFZR_001274(self):
        """
        This test takes a blast result and return ltg_tax_id, ltg_rank at given identity
        """
        #
        blast_out_df = f04_import_blast_output_into_df(self.blast_MFZR_001274_tsv)
        #
        identity = 80
        blast_out_df.loc[blast_out_df.identity >= self.identity_threshold]
        blast_result_subset_df = blast_out_df.loc[blast_out_df.identity >= identity, ['target_id', 'target_tax_id']]
        tax_lineage_df = f05_blast_result_subset(blast_result_subset_df, self.taxonomy_db_df)
        ltg_tax_id, ltg_rank = f06_select_ltg(tax_lineage_df, identity=identity, identity_threshold=self.identity_threshold)
        self.assertTrue(ltg_tax_id==1344033)
        self.assertTrue(ltg_rank=='species')

    def test_f07_blast_result_to_ltg(self):
        blast_result_to_ltg = """13	HF558455.1	86.471	170	23	0	1	170	9324	9493	1.91e-48	204
13	MF694646.1	83.529	170	28	0	1	170	11574	11743	6.25e-42	181
13	KY911088.1	80.791	177	34	0	1	177	31380	31556	1.38e-37	167
13	KY911087.1	80.791	177	34	0	1	177	31501	31677	1.38e-37	167
13	JQ007431.1	80.791	177	34	0	1	177	32	208	1.38e-37	167
13	JQ007430.1	80.791	177	34	0	1	177	32	208	1.38e-37	167
13	JQ007428.1	80.791	177	34	0	1	177	32	208	1.38e-37	167
13	KY911089.1	80.226	177	35	0	1	177	32943	33119	5.86e-36	162
13	JQ007429.1	80.226	177	35	0	1	177	32	208	5.86e-36	162
13	LK052976.1	78.212	179	39	0	2	180	8577	8755	1.29e-31	148
13	DQ157700.1	78.035	173	38	0	1	173	30341	30169	5.49e-30	141
13	KC628747.1	77.222	180	41	0	1	180	13088	12909	1.91e-29	141
13	HG529787.1	77.778	171	38	0	1	171	74574	74404	6.68e-29	138
13	NC_039443.1	77.457	173	39	0	2	174	43	215	2.33e-28	137
13	MH725800.1	77.457	173	39	0	2	174	43	215	2.33e-28	137
13	MG783568.1	77.011	174	40	0	2	175	210	383	8.14e-28	134
13	AP017979.1	77.711	166	37	0	3	168	8656	8821	2.84e-27	133
13	CP010939.1	77.059	170	39	0	1	170	77066	76897	9.92e-27	132
13	KY911092.1	76.744	172	40	0	2	173	30226	30397	9.92e-27	131
13	XM_018135786.1	76.744	172	40	0	2	173	46	217	9.92e-27	131
13	MH725798.1	76.437	174	41	0	2	175	282	455	3.46e-26	130
13	JX499144.1	76.301	173	41	0	3	175	18809	18981	1.21e-25	128
13	KY911091.1	75.706	177	43	0	2	178	14673	14849	4.22e-25	126
13	XM_012330695.1	79.577	142	29	0	32	173	295	436	4.22e-25	126
13	NC_039442.1	75.419	179	44	0	2	180	4268	4446	4.22e-25	125
13	MH725794.1	75.419	179	44	0	2	180	4268	4446	4.22e-25	125
13	KX810692.1	81.746	126	23	0	56	181	17	142	1.47e-24	124
13	KX810689.1	81.746	126	23	0	56	181	17	142	1.47e-24	124
13	KX810686.1	83.621	116	19	0	56	171	17	132	1.47e-24	124
13	AP017981.1	76.506	166	39	0	3	168	29783	29948	1.47e-24	124
13	AP017980.1	76.506	166	39	0	3	168	35698	35863	1.47e-24	124
13	NC_039441.1	75.723	173	42	0	2	174	12053	12225	1.47e-24	123
13	MH725793.1	75.723	173	42	0	2	174	12053	12225	1.47e-24	123
13	FQ311469.1	75.740	169	41	0	2	170	57966	57798	1.79e-23	121
13	NC_039439.1	75.145	173	43	0	2	174	79	251	6.26e-23	119
13	MH725791.1	75.145	173	43	0	2	174	79	251	6.26e-23	119
13	KY245891.1	76.074	163	39	0	3	165	47251	47413	6.26e-23	119
13	XM_025501004.1	75.294	170	42	0	2	171	46	215	6.26e-23	118
13	LC385608.1	76.506	166	37	2	11	175	28217	28381	2.18e-22	117
13	KY911097.1	76.437	174	35	4	1	171	34051	34221	2.18e-22	117
13	JX985789.1	76.506	166	37	2	11	175	52	216	2.18e-22	117
13	NC_039440.1	74.713	174	44	0	2	175	140579	140406	2.18e-22	116
13	MH725792.1	74.713	174	44	0	2	175	140579	140406	2.18e-22	116
13	JX499145.1	75.148	169	42	0	3	171	20046	20214	2.18e-22	116
13	GU133622.1	75.148	169	42	0	3	171	4888	4720	2.18e-22	116
13	KY911085.1	74.854	171	43	0	1	171	11396	11566	7.63e-22	115
13	KY911084.1	74.854	171	43	0	1	171	11560	11730	7.63e-22	115
13	KY911083.1	74.854	171	43	0	1	171	11537	11707	7.63e-22	115
13	KY911082.1	74.854	171	43	0	1	171	11342	11512	7.63e-22	115
13	KY911081.1	74.854	171	43	0	1	171	11311	11481	7.63e-22	115"""
        df = pandas.DataFrame([x.split(';') for x in data.split('\n')])