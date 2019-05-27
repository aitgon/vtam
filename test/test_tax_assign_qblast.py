# -*- coding: utf-8 -*-
import inspect
import os
import sqlite3
import urllib

import pandas
from unittest import TestCase

from wopmetabarcoding.utils.PathFinder import PathFinder
from wopmetabarcoding.utils.constants import VTAM_DATA_DIR, public_data_dir
from wopmetabarcoding.utils.logger import logger
from wopmetabarcoding.wrapper.TaxAssignUtilities import f01_taxonomy_sqlite_to_df, f06_select_ltg, \
    f04_import_qblast_output_into_df, f05_blast_result_subset


class TestTaxAssign(TestCase):

    def setUp(self):
        #####################################
        #
        # Download taxonomy.sqlite
        #
        #####################################

        PathFinder.mkdir_p(VTAM_DATA_DIR)
        file_remote = os.path.join(public_data_dir, "taxonomy.sqlite")
        taxonomy_sqlite_path = os.path.join(os.environ['DIR_DATA_NON_GIT'], 'taxonomy.sqlite')
        # taxonomy_sqlite_path = os.path.join(VTAM_DATA_DIR, os.path.basename(file_remote))
        if not os.path.isfile(taxonomy_sqlite_path):
            logger.debug(
                "file: {}; line: {}; Downloading taxonomy.sqlite".format(__file__, inspect.currentframe().f_lineno))
            urllib.request.urlretrieve(file_remote, taxonomy_sqlite_path)
        #
        self.taxonomy_db_df = f01_taxonomy_sqlite_to_df(taxonomy_sqlite_path)
        #
        self.identity_threshold = 97
        self.min_number_of_taxa = 3
        self.include_prop = 90
        #
        self.__testdir_path = os.path.join(PathFinder.get_module_test_path())
        self.qblast_MFZR_002737_tsv = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "test_files", "qblast_MFZR_002737.tsv")
        self.qblast_MFZR_001274_tsv = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "test_files", "qblast_MFZR_001274.tsv")
        self.v1_fasta = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "MFZR_001274.fasta")
        self.qblast_MFZR_002737_tsv = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "test_files","qblast_MFZR_002737.tsv")
        self.tax_lineage_variant13_tsv = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "test_files","tax_lineage_variant13.tsv")
        # self.tax_lineage_tsv = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path,
        #                                               "test_files", "tax_lineage.tsv")

    #
    # Commented because too slow
    # Uncomment to test qqblast
    # def test_f06_2_run_qblast(self):
    #     #
    #     # Run and read qblast result
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

    def test_f06_3_annotate_qblast_output_with_tax_id(self):
        #
        qblast_MFZR_001274_tsv = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "test_files",
                                                  "qblast_MFZR_001274.tsv")
        nucl_gb_accession2taxid_MFZR_00001274_sqlite = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "test_files",
                                                  "nucl_gb_accession2taxid_MFZR_001274.sqlite")
        #
        # Run and read qblast result
        qblast_result_df = pandas.read_csv(qblast_MFZR_001274_tsv, sep="\t", skiprows=13, usecols=[0, 1, 2],
                             header=None, names=['variant_id', 'gb_accession', 'identity'])
        qblast_result_df = qblast_result_df.loc[~pandas.isnull(qblast_result_df).any(axis=1)]
        # import pdb; pdb.set_trace()
        con = sqlite3.connect(nucl_gb_accession2taxid_MFZR_00001274_sqlite)
        sql = """SELECT gb_accession, tax_id FROM nucl_gb_accession2taxid WHERE gb_accession IN {}""".format(tuple(qblast_result_df.gb_accession.tolist()))
        gb_accession_to_tax_id_df = pandas.read_sql(sql=sql, con=con)
        con.close()
        #
        qblast_result_tax_id_df = qblast_result_df.merge(gb_accession_to_tax_id_df, on='gb_accession')


    def test_f03_import_qblast_output_into_df(self):
        """This test assess whether the qblast has been well imported into a df"""
        #
        qblast_out_df = f04_import_qblast_output_into_df(self.qblast_MFZR_002737_tsv)
        #
        self.assertTrue(qblast_out_df.to_dict('list') == {'target_id': [1049499563, 1049496963, 1049491687, 1049490545],
                                         'identity': [100.0, 100.0, 99.429, 99.429],
                                         'target_tax_id': [761875, 761875, 761875, 761875]})

    def test_99_full_tax_assign_after_qblast_MFZR_002737(self):
        """
        This test takes a qblast result and return ltg_tax_id, ltg_rank at given identity
        """
        #
        qblast_out_df = f04_import_qblast_output_into_df(self.qblast_MFZR_002737_tsv)
        #
        identity = 100
        #
        qblast_out_df.loc[qblast_out_df.identity >= self.identity_threshold]
        qblast_result_subset_df = qblast_out_df.loc[qblast_out_df.identity >= identity, ['target_id', 'target_tax_id']]
        tax_lineage_df = f05_blast_result_subset(qblast_result_subset_df, self.taxonomy_db_df)
        ltg_tax_id, ltg_rank = f06_select_ltg(tax_lineage_df=tax_lineage_df, identity=identity,
                            identity_threshold=self.identity_threshold, include_prop=self.include_prop,
                                              min_number_of_taxa=self.min_number_of_taxa)
        self.assertTrue(ltg_tax_id==761875)
        self.assertTrue(ltg_rank=='species')

    def test_99_full_tax_assign_after_qblast_MFZR_001274(self):
        """
        This test takes a qblast result and return ltg_tax_id, ltg_rank at given identity
        """
        #
        qblast_out_df = f04_import_qblast_output_into_df(self.qblast_MFZR_001274_tsv)
        #
        identity = 80
        qblast_out_df.loc[qblast_out_df.identity >= self.identity_threshold]
        qblast_result_subset_df = qblast_out_df.loc[qblast_out_df.identity >= identity, ['target_id', 'target_tax_id']]
        tax_lineage_df = f05_blast_result_subset(qblast_result_subset_df, self.taxonomy_db_df)
        ltg_tax_id, ltg_rank = f06_select_ltg(tax_lineage_df=tax_lineage_df, identity=identity,
                            identity_threshold=self.identity_threshold, include_prop=self.include_prop,
                                              min_number_of_taxa=self.min_number_of_taxa)
        self.assertTrue(ltg_tax_id==1344033)
        self.assertTrue(ltg_rank=='species')
