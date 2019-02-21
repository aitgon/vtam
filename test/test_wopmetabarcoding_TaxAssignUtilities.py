import os
from unittest import TestCase

import pandas
from numpy import nan

from wopmetabarcoding.utils.PathFinder import PathFinder
from wopmetabarcoding.wrapper.TaxAssignUtilities import create_phylogenetic_line_df, f_majoritytaxid2percentage, \
    f_taxid2taxname, f_taxlineage_to_ltg


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


    # def test_create_phylogenetic_line_df(self):
    #     test_outdir = os.path.join(self.__testdir_path, "output", "TaxAssignUtilities")
    #     PathFinder.mkdir_p(test_outdir)
    #     tax_seq_id_list_txt = os.path.join(PathFinder.get_module_test_path(), "input", "TaxAssignUtilities", "tax_seq_id_list.txt")
    #     tax_lineage_df_pkl = os.path.join(PathFinder.get_module_test_path(), "input", "TaxAssignUtilities", "tax_lineage_df.pkl")
    #     #
    #     # Input
    #     with open(tax_seq_id_list_txt, 'r') as fin:
    #         tax_seq_id_list = [l.strip() for l in fin]
    #     #
    #     # Process
    #     tax_lineage_df = create_phylogenetic_line_df(tax_seq_id_list, self.TAX_ASSIGN_SQLITE)
    #     #
    #     # Output
    #     tax_lineage_df_bak = pandas.read_pickle(tax_lineage_df_pkl)
    #     #
    #     tax_lineage_df_bak.sort_values(by=tax_lineage_df_bak.columns.tolist(), inplace=True)
    #     tax_lineage_df.sort_values(by=tax_lineage_df.columns.tolist(), inplace=True)
    #     tax_lineage_df = tax_lineage_df.apply(pandas.to_numeric)
    #     tax_lineage_df_bak = tax_lineage_df_bak.apply(pandas.to_numeric)
    #     tax_lineage_df.equals(tax_lineage_df_bak)
    #     self.assertTrue(tax_lineage_df.equals(tax_lineage_df_bak))
    #
    #
    # def test_dataframe2ltgdefinition(self):
    #     test_outdir = os.path.join(self.__testdir_path, "output", "TaxAssignUtilities")
    #     PathFinder.mkdir_p(test_outdir)
    #     tax_count_perc_df_pkl = os.path.join(PathFinder.get_module_test_path(), "input", "TaxAssignUtilities", "tax_count_perc_df.pkl")
    #     #
    #     # Input
    #     tax_lineage_df_pkl = os.path.join(PathFinder.get_module_test_path(), "input", "TaxAssignUtilities", "tax_lineage_df.pkl")
    #     tax_lineage_df = pandas.read_pickle(tax_lineage_df_pkl)
    #     #
    #     tax_count_perc_bak_df = pandas.read_pickle(tax_count_perc_df_pkl)
    #     #
    #     # Run
    #     tax_count_perc_df = f_majoritytaxid2percentage(tax_lineage_df)
    #     self.assertTrue(tax_count_perc_df.equals(tax_count_perc_bak_df))
    #
    #
    # def test_f_taxid2taxname(self):
    #     tax_id_list = [603345, 325869, 0]
    #     taxid2taxname_dic = f_taxid2taxname(tax_id_list, "/home/gonzalez/Data/2017_meglecz_metabarcoding/data/wopmetabarcoding_data/tax_assign.sqlite")
    #     self.assertTrue(taxid2taxname_dic == {603345: 'Zebrias synapturoides', 325869: 'Xylophanes tersa', 0: nan})
    #
    #
    # def test_f_taxlineage_to_ltg(self):
    #     # test_outdir = os.path.join(self.__testdir_path, "output", "TaxAssignUtilities")
    #     # PathFinder.mkdir_p(test_outdir)
    #     tax_count_perc_df_pkl = os.path.join(PathFinder.get_module_test_path(), "input", "TaxAssignUtilities", "tax_count_perc_df.pkl")
    #     #
    #     # Input
    #     tax_lineage_df_pkl = os.path.join(PathFinder.get_module_test_path(), "input", "TaxAssignUtilities", "f_taxlineage_to_ltg", "tax_lineage_df.pkl")
    #     tax_lineage_df = pandas.read_pickle(tax_lineage_df_pkl)
    #     #
    #     max_tax_resolution_id = 15
    #     self.assertTrue(f_taxlineage_to_ltg(tax_lineage_df, max_tax_resolution_id), 325869)
