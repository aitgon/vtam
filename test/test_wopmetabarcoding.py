import filecmp
import os
import pickle
import shutil
import sqlite3
from unittest import TestCase
import subprocess

import pandas
from jinja2 import Template
from wopmars import WopMars
from wopmars.framework.parsing.Parser import ToolWrapper

from wopmetabarcoding.utils.PathFinder import PathFinder
from wopmetabarcoding.wrapper.FilterUtilities import Variant2Sample2Replicate2Count
from wopmetabarcoding.utils.constants import wopmetabarcoding_filter_test_data



class TestWopMetabarcoding(TestCase):
    def setUp(self):
        self.__testdir_path = os.path.join(PathFinder.get_module_test_path())
        self.__db_path = os.path.join(PathFinder.get_module_test_path(), self.__testdir_path, "db.sqlite")
        self.__db_url = "sqlite:///" + self.__db_path
        self.__wopfile_template_path = PathFinder.get_wopfile_template_path()
        self.__wopfile_template_str = ""
        with open(self.__wopfile_template_path, 'r') as fin:
            self.__wopfile_template_str = fin.read()

    def test_01_sample_information(self):
        test_outdir = os.path.join(self.__testdir_path, "output", "01sample_information")
        PathFinder.mkdir(test_outdir)
        sample_info_csv = os.path.join(PathFinder.get_module_test_path(), "input", "01sample_information", "1_sample_info.csv")
        db_path = os.path.join(test_outdir, "db.sqlite")
        db_url = "sqlite:///" + db_path
        #
        # Create wopfile
        template = Template(self.__wopfile_template_str)
        wopfile_str = template.render(SAMPLE_INFORMATION_CSV=sample_info_csv)
        wopfile_path = os.path.join(test_outdir, "Wopfile.yml")
        with open(wopfile_path, 'w') as fout:
            self.__wopfile_template_str = fout.write(wopfile_str)
        #
        cmd_line = ["wopmars", "-w", wopfile_path, "-D", db_url, "-d", self.__testdir_path, "-F", "-t", "SampleInformation"]
        with self.assertRaises(SystemExit) as se:
            WopMars().run(cmd_line)
        self.assertEqual(se.exception.code, 0)
        #
        output = (1, 'ZFZR', 'AGATATTGGAACWTTATATTTTATTTTTGG', 'WACTAATCAATTWCCAAATCCTCC', 'cgatcgtcatcacg', 'ctcgatgatcacg', 1, 'prerun', '14Mon01', 'repl1')
        con = sqlite3.connect(self.__db_path)
        cur = con.cursor()
        cur.execute("SELECT * from SampleInformation;")
        self.assertTrue(list(cur.fetchone()) == list(output))
        #
        shutil.rmtree(test_outdir)

    def test_02_sample_information_error(self):
        test_outdir = os.path.join(self.__testdir_path, "output", "01sample_information")
        PathFinder.mkdir(test_outdir)
        sample_info_csv = os.path.join(PathFinder.get_module_test_path(), "input", "01sample_information", "1_sample_info_error.csv")
        db_path = os.path.join(test_outdir, "db.sqlite")
        db_url = "sqlite:///" + db_path
        #
        # Create wopfile
        template = Template(self.__wopfile_template_str)
        wopfile_str = template.render(SAMPLE_INFORMATION_CSV=sample_info_csv)
        wopfile_path = os.path.join(test_outdir, "Wopfile.yml")
        with open(wopfile_path, 'w') as fout:
            self.__wopfile_template_str = fout.write(wopfile_str)
        #
        cmd_line = ["wopmars", "-w", wopfile_path, "-D", db_url, "-d", self.__testdir_path, "-F", "-t", "SampleInformation"]
        with self.assertRaises(SystemExit) as se:
            WopMars().run(cmd_line)
        self.assertEqual(se.exception.code, 1)
        #
        shutil.rmtree(test_outdir)

    def test_03_sort_reads(self):
        test_outdir = os.path.join(self.__testdir_path, "output", "02sort_reads")
        PathFinder.mkdir(test_outdir)
        sample_info_csv = os.path.join(PathFinder.get_module_test_path(), "input", "01sample_information", "1_sample_info.csv")
        db_path = os.path.join(test_outdir, "db.sqlite")
        db_url = "sqlite:///" + db_path
        #
        # Create wopfile
        template = Template(self.__wopfile_template_str)
        wopfile_str = template.render(SAMPLE_INFORMATION_CSV=sample_info_csv)
        wopfile_path = os.path.join(test_outdir, "Wopfile.yml")
        with open(wopfile_path, 'w') as fout:
            self.__wopfile_template_str = fout.write(wopfile_str)
        #
        cmd_line = ["wopmars", "-w", wopfile_path, "-D", db_url, "-d", self.__testdir_path, "-F", "-t", "SortReads"]
        p = subprocess.Popen(cmd_line)
        p.wait()
        #
        marker_id = 1
        sequence = 'AGCCTGAGCTGGAATAGTAGGTACTTCCCTTAGTATACTTATTCGAGCCGAATTAGGACACCCAGGCTCTCTAATTGGAGACGACCAAATTTATAATGTAATTGTTACTGCTCATGCTTTTGTAATAATTTTTTTTATAGTTATGCCAATTATAATT'
        con = sqlite3.connect(db_path)
        cur = con.cursor()
        sql = "SELECT marker_id from Variant where sequence='%s';"%(sequence)
        cur.execute(sql)
        self.assertTrue(list(cur.fetchone()) == [1])
        #
        shutil.rmtree(test_outdir)

    def test_04_store_index_below_lfn1_per_replicate(self):
        test_outdir = os.path.join(self.__testdir_path, "output", "03filter")
        PathFinder.mkdir(test_outdir)
        variant2sample2replicate2count_df_pkl_path = os.path.join(PathFinder.get_module_test_path(), "input", "03filter", "variant2sample2replicate2count_df.pkl")
        #
        # Input
        variant2sample2replicate2count_df = pandas.read_pickle(variant2sample2replicate2count_df_pkl_path)
        lfn1_per_replicate_threshold = 0.25
        #
        # Output
        variant2sample2replicate2count = Variant2Sample2Replicate2Count(variant2sample2replicate2count_df)
        variant2sample2replicate2count.store_index_below_lfn1_per_replicate(lfn1_per_replicate_threshold)
        indices_to_drop = variant2sample2replicate2count.indices_to_drop
        indices_to_drop_bak = [76, 88, 117, 119]
        self.assertTrue(indices_to_drop == indices_to_drop_bak)
        #
        shutil.rmtree(test_outdir)

    def test_05_store_index_below_lfn2_per_variant(self):
        test_outdir = os.path.join(self.__testdir_path, "output", "03filter")
        PathFinder.mkdir(test_outdir)
        variant2sample2replicate2count_df_pkl_path = os.path.join(PathFinder.get_module_test_path(), "input", "03filter", "variant2sample2replicate2count_df.pkl")
        #
        # Input
        variant2sample2replicate2count_df = pandas.read_pickle(variant2sample2replicate2count_df_pkl_path)
        lfn2_per_variant_threshold = 0.025
        #
        # Output
        variant2sample2replicate2count = Variant2Sample2Replicate2Count(variant2sample2replicate2count_df)
        variant2sample2replicate2count.store_index_below_lfn2_per_variant(lfn2_per_variant_threshold)
        indices_to_drop = variant2sample2replicate2count.indices_to_drop
        # print(indices_to_drop)
        indices_to_drop_bak = [36, 46, 53]
        self.assertTrue(indices_to_drop == indices_to_drop_bak)
        #
        shutil.rmtree(test_outdir)

    def test_06_store_index_below_lfn2_per_replicate_series(self):
        test_outdir = os.path.join(self.__testdir_path, "output", "03filter")
        PathFinder.mkdir(test_outdir)
        variant2sample2replicate2count_df_pkl_path = os.path.join(PathFinder.get_module_test_path(), "input", "03filter", "variant2sample2replicate2count_df.pkl")
        #
        # Input
        variant2sample2replicate2count_df = pandas.read_pickle(variant2sample2replicate2count_df_pkl_path)
        lfn2_per_variant_threshold = 0.0075
        #
        # Output
        variant2sample2replicate2count = Variant2Sample2Replicate2Count(variant2sample2replicate2count_df)
        variant2sample2replicate2count.store_index_below_lfn2_per_replicate_series(lfn2_per_variant_threshold)
        indices_to_drop = variant2sample2replicate2count.indices_to_drop
        # print(indices_to_drop)
        indices_to_drop_bak = [27, 36, 46, 53, 88, 92, 104, 122, 209]
        self.assertTrue(indices_to_drop == indices_to_drop_bak)
        #
        shutil.rmtree(test_outdir)

    def test_07_store_index_below_lfn3_read_count(self):
        test_outdir = os.path.join(self.__testdir_path, "output", "03filter")
        PathFinder.mkdir(test_outdir)
        variant2sample2replicate2count_df_pkl_path = os.path.join(PathFinder.get_module_test_path(), "input", "03filter", "variant2sample2replicate2count_df.pkl")
        #
        # Input
        variant2sample2replicate2count_df = pandas.read_pickle(variant2sample2replicate2count_df_pkl_path)
        lfn_read_count_threshold = 3
        #
        # Output
        variant2sample2replicate2count = Variant2Sample2Replicate2Count(variant2sample2replicate2count_df)
        variant2sample2replicate2count.store_index_below_lfn3_read_count(lfn_read_count_threshold)
        indices_to_drop = variant2sample2replicate2count.indices_to_drop
        # print(indices_to_drop)
        indices_to_drop_bak = [27, 36, 46, 53, 88, 92, 104, 122, 209]
        self.assertTrue(indices_to_drop == indices_to_drop_bak)
        #
        shutil.rmtree(test_outdir)

    # def test_class_variant2sample2replicate2count(self):
    #     #
    #     # Input
    #     db_sqlite = os.path.join(wopmetabarcoding_filter_test_data, "db.sqlite")
    #     ZFZR_sample_count_tsv = os.path.join(wopmetabarcoding_filter_test_data, "ZFZR_sample_count.tsv")
    #     variant2sample2replicate2count_df = pandas.read_csv(ZFZR_sample_count_tsv, sep='\t')
    #     variant2sample2replicate2count = Variant2Sample2Replicate2Count(variant2sample2replicate2count_df)
    #     lfn_per_replicate_threshold = 0.025
    #     failed_indices = [309, 353, 55, 237, 377, 263, 354, 408, 145, 178, 240, 287, 319, 339]
    #     self.assertTrue(variant2sample2replicate2count.store_index_below_lfn1_per_replicate(lfn_per_replicate_threshold) == failed_indices)


    # def test_03_chimera(self):
    #     df_pkl = os.path.join(PathFinder.get_module_test_path(), "input/filter/store_index_identified_as_chimera", "df.pkl")
    #     replicate_obj_list_pkl = os.path.join(PathFinder.get_module_test_path(), "input/filter/store_index_identified_as_chimera", "replicate_obj_list.pkl")
    #     with open(df_pkl, 'rb') as fin:
    #         df = pickle.load(fin)
    #     with open(replicate_obj_list_pkl, 'rb') as fin:
    #         replicate_obj_list = pickle.load(fin)
    #     marker_id = 1
    #     chimera_by = "sample_replicate"
    #     store_index_identified_as_chimera(replicate_obj_list, df, marker_id, chimera_by)

