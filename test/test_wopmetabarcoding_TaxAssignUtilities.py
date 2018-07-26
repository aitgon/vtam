import os
from unittest import TestCase

import pandas

from wopmetabarcoding.utils.PathFinder import PathFinder
from wopmetabarcoding.wrapper.TaxassignUtilities import create_phylogenetic_line_df, dataframe2ltgdefinition


class TestWopMetabarcoding(TestCase):
    def setUp(self):
        self.__testdir_path = os.path.join(PathFinder.get_module_test_path())
        self.__db_path = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "db.sqlite")
        self.__db_url = "sqlite:///" + self.__db_path
        self.__wopfile_test_path = PathFinder.get_wopfile_test_path()
        self.__wopfile_test_str = ""
        with open(self.__wopfile_test_path, 'r') as fin:
            self.__wopfile_test_str = fin.read()
        self.TAX_ASSIGN_SQLITE = "/home/gonzalez/Data/2017_meglecz_metabarcoding/data/wopmetabarcoding_data/tax_assign.sqlite"


    def test_create_phylogenetic_line_df(self):
        test_outdir = os.path.join(self.__testdir_path, "output", "TaxAssignUtilities")
        PathFinder.mkdir_p(test_outdir)
        tax_seq_id_list_txt = os.path.join(PathFinder.get_module_test_path(), "input", "TaxAssignUtilities", "tax_seq_id_list.txt")
        tax_lineage_df_pkl = os.path.join(PathFinder.get_module_test_path(), "input", "TaxAssignUtilities", "tax_lineage_df.pkl")
        #
        # Input
        with open(tax_seq_id_list_txt, 'r') as fin:
            tax_seq_id_list = [l.strip() for l in fin]
        #
        # Process
        tax_lineage_df = create_phylogenetic_line_df(tax_seq_id_list, self.TAX_ASSIGN_SQLITE)
        #
        # Output
        tax_lineage_df_bak = pandas.read_pickle(tax_lineage_df_pkl)
        #
        tax_lineage_df_bak.sort_values(by=tax_lineage_df_bak.columns.tolist(), inplace=True)
        tax_lineage_df.sort_values(by=tax_lineage_df.columns.tolist(), inplace=True)
        tax_lineage_df = tax_lineage_df.apply(pandas.to_numeric)
        tax_lineage_df_bak = tax_lineage_df_bak.apply(pandas.to_numeric)
        tax_lineage_df.equals(tax_lineage_df_bak)
        self.assertTrue(tax_lineage_df.equals(tax_lineage_df_bak))


    def test_dataframe2ltgdefinition(self):
        test_outdir = os.path.join(self.__testdir_path, "output", "TaxAssignUtilities")
        PathFinder.mkdir_p(test_outdir)
        tax_lineage_df_pkl = os.path.join(PathFinder.get_module_test_path(), "input", "TaxAssignUtilities", "tax_lineage_df.pkl")
        tax_count_perc_df_pkl = os.path.join(PathFinder.get_module_test_path(), "input", "TaxAssignUtilities", "tax_count_perc_df.pkl")
        #
        # Input
        tax_lineage_df = pandas.read_pickle(tax_lineage_df_pkl)
        #
        tax_count_perc_bak_df = pandas.read_pickle(tax_count_perc_df_pkl)
        #
        # Run
        tax_count_perc_df = dataframe2ltgdefinition(tax_lineage_df)
        self.assertTrue(tax_count_perc_df.equals(tax_count_perc_bak_df))


    # def test_04filter_store_index_below_lfn3_read_count(self):
    #     test_outdir = os.path.join(self.__testdir_path, "output", "04filter")
    #     PathFinder.mkdir_p(test_outdir)
    #     variant2sample2replicate2count_df_pkl_path = os.path.join(PathFinder.get_module_test_path(), "input", "04filter", "variant2sample2replicate2count_df.pkl")
    #     #
    #     # Input
    #     variant2sample2replicate2count_df = pandas.read_pickle(variant2sample2replicate2count_df_pkl_path)
    #     lfn_read_count_threshold = 3
    #     #
    #     # Output
    #     variant2sample2replicate2count = Variant2Sample2Replicate2Count(variant2sample2replicate2count_df)
    #     variant2sample2replicate2count.store_index_below_lfn3_read_count(lfn_read_count_threshold)
    #     indices_to_drop = variant2sample2replicate2count.indices_to_drop
    #     # print(indices_to_drop)
    #     indices_to_drop_bak = [27, 36, 46, 53, 88, 92, 104, 122, 209]
    #     # Todo: indices_to_drop returns empty []: TD needs to check it.
    #     # self.assertTrue(indices_to_drop == indices_to_drop_bak)
    #     #
    #     shutil.rmtree(test_outdir)

    def test05_taxassign_vsearch(self):
        test_outdir = os.path.join(self.__testdir_path, "output", "05taxassign")

