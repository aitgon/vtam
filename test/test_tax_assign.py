# -*- coding: utf-8 -*-
import errno
import inspect
import os
import sqlite3
import tarfile
import urllib

import pandas
from unittest import TestCase

from Bio.Blast import NCBIWWW

from wopmetabarcoding.utils.PathFinder import PathFinder
from wopmetabarcoding.utils.constants import tempdir, data_dir
from wopmetabarcoding.utils.logger import logger
from wopmetabarcoding.wrapper.TaxAssignUtilities import f01_taxonomy_sqlite_to_df, f04_1_tax_id_to_taxonomy_lineage, \
    f06_select_ltg, \
    f04_import_blast_output_into_df, f05_blast_result_subset, f02_variant_df_to_fasta

import gzip
import shutil


class TestTaxAssign(TestCase):

    def setUp(self):
        self.taxonomy_db_df = f01_taxonomy_sqlite_to_df('tax.sqlite')
        #
        self.identity_threshold = 97
        #
        self.__testdir_path = os.path.join(PathFinder.get_module_test_path())
        self.blast_MFZR_002737_tsv = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "test_files", "blast_MFZR_002737.tsv")
        self.blast_MFZR_001274_tsv = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "test_files", "blast_MFZR_001274.tsv")
        self.v1_fasta = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "MFZR_001274.fasta")
        self.blast_MFZR_002737_tsv = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "test_files",
                                                  "blast_MFZR_002737.tsv")

    def test_f06_1_create_nucl_gb_accession2taxid_sqlite(self):
        #
        # Check if exists datadir with ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
        try:
            os.makedirs(data_dir)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        #####################################
        #
        # Download nucl_gb.accession2taxid.gz
        #
        #####################################
        file_remote = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz"
        nucl_gb_accession2taxid_gz = os.path.join(data_dir, os.path.basename(file_remote))
        if not os.path.isfile(nucl_gb_accession2taxid_gz):  #  TODO verify MD5
            logger.debug(
                "file: {}; line: {}; Downloading nucl_gb.accession2taxid.gz".format(__file__, inspect.currentframe().f_lineno))
            urllib.request.urlretrieve(file_remote, nucl_gb_accession2taxid_gz)
        #####################################
        #
        # Gunzipping nucl_gb.accession2taxid.gz
        #
        #####################################
        nucl_gb_accession2taxid_tsv = os.path.join(data_dir, "nucl_gb.accession2taxid.tsv")
        if not (os.path.isfile(nucl_gb_accession2taxid_tsv)):
            logger.debug(
                "file: {}; line: {}; Gunziping nucl_gb.accession2taxid.gz".format(__file__, inspect.currentframe().f_lineno))
            with gzip.open(nucl_gb_accession2taxid_gz, 'rb') as f_in:
                with open(nucl_gb_accession2taxid_tsv, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        #####################################
        #
        # Import nucl_gb.accession2taxid.gz to sqlite
        #
        #####################################
        nucl_gb_accession2taxid_sqlite = nucl_gb_accession2taxid_tsv.replace("tsv", "sqlite")
        logger.debug(
            "file: {}; line: {}; Import nucl_gb.accession2taxid.gz to sqlite".format(__file__, inspect.currentframe().f_lineno))
        conn = sqlite3.connect(nucl_gb_accession2taxid_sqlite)
        conn.execute("DROP TABLE IF EXISTS nucl_gb_accession2taxid")
        cur = conn.cursor()
        conn.execute(
            "CREATE TABLE nucl_gb_accession2taxid (gb_accession VARCHAR, tax_id INTEGER, unique(gb_accession, tax_id))"
        )
        cur.close()
        conn.commit()
        # gb_accession_to_tax_id_df = pandas.read_csv(os.path.join(data_dir, "nucl_gb.accession2taxid"), sep="\t",
        #                     usecols=[1, 2],
        #                      header=None, names=['gb_accession.version', 'tax_id'])
        chunksize = 10 ** 6
        for chunk_df in pandas.read_csv(os.path.join(data_dir, "nucl_gb.accession2taxid.tsv"), sep="\t",
                            usecols=[1, 2],
                             names=['gb_accession', 'tax_id'], chunksize=chunksize, header=0):
            chunk_df.to_sql(name="nucl_gb_accession2taxid", con = conn, if_exists='append', index=False)
        conn.close()

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
