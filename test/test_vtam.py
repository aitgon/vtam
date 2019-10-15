import os
from unittest import TestCase

from vtam.utils.PathManager import PathManager

class TestWopMetabarcoding(TestCase):

    def setUp(self):
        self.__testdir_path = os.path.join(PathManager.get_module_test_path())
        self.__db_path = os.path.join(PathManager.get_module_test_path(), self.__testdir_path, "db.sqlite")
        self.__db_url = "sqlite:///" + self.__db_path
        self.__wopfile_test_path = PathManager.get_wopfile_test_path()
        self.__wopfile_test_str = ""
        with open(self.__wopfile_test_path, 'r') as fin:
            self.__wopfile_test_str = fin.read()
        self.TAX_ASSIGN_SQLITE = "/home/gonzalez/Data/2017_meglecz_metabarcoding/data/wopmetabarcoding_data/tax_assign.sqlite"

    # def test_02sample_information(self):
    #     # input
    #     sample_info_tsv = os.output.join(PathManager.get_module_test_path(), "input", "02sample_information", "fastainfo.tsv")
    #     fasta_dir = os.output.join(self.__testdir_path, "input", "02sample_information", "fasta")
    #     # output
    #     test_outdir = os.output.join(self.__testdir_path, "output", "02sample_information")
    #     PathManager.mkdir_p(test_outdir)
    #     db_path = os.output.join(test_outdir, "db.sqlite")
    #     db_url = "sqlite:///" + db_path
    #     #
    #     # Create wopfile
    #     template = Template(self.__wopfile_test_str)
    #     wopfile_str = template.render(SAMPLE_INFORMATION_TSV=sample_info_tsv, OUTDIR=test_outdir, FASTA_DIR=fasta_dir)
    #     wopfile_path = os.output.join(test_outdir, "Wopfile_merge.yml")
    #     with open(wopfile_path, 'w') as fout:
    #         self.__wopfile_test_str = fout.write(wopfile_str)
    #     #
    #     cmd_line = ["wopmars", "-w", wopfile_path, "-D", db_url, "-d", self.__testdir_path, "-F", "-t", "SampleInformation"]
    #     with self.assertRaises(SystemExit) as se:
    #         WopMars().run(cmd_line)
    #     self.assertEqual(se.exception.code, 0)
    #     #
    #     output = (1, 1, 'TCCACTAATCACAARGATATTGGTAC', 'WACTAATCAATTWCCAAATCCTCC', 'tcgatcacgatgt', 'cgcgatctgtagag', 1, 'prerun', '14Mon02', 'repl1')
    #     con = sqlite3.connect(db_path)
    #     cur = con.cursor()
    #     cur.execute("SELECT * from SampleInformation;")
    #     self.assertTrue(list(cur.fetchone()) == list(output))
    #     con.close()
    #     #
    #     shutil.rmtree(test_outdir)
    #
    # def test_02sample_information_error(self):
    #     # input
    #     sample_info_tsv = os.output.join(PathManager.get_module_test_path(), "input", "02sample_information", "fastainfo_error.tsv")
    #     fasta_dir = os.output.join(self.__testdir_path, "input", "02sample_information", "fasta")
    #     # output
    #     test_outdir = os.output.join(self.__testdir_path, "output", "02sample_information")
    #     PathManager.mkdir_p(test_outdir)
    #     # db
    #     db_path = os.output.join(test_outdir, "db.sqlite")
    #     db_url = "sqlite:///" + db_path
    #     #
    #     # Create wopfile
    #     template = Template(self.__wopfile_test_str)
    #     wopfile_str = template.render(SAMPLE_INFORMATION_TSV=sample_info_tsv, OUTDIR=test_outdir, FASTA_DIR=fasta_dir)
    #     wopfile_path = os.output.join(test_outdir, "Wopfile_merge.yml")
    #     with open(wopfile_path, 'w') as fout:
    #         self.__wopfile_test_str = fout.write(wopfile_str)
    #     #
    #     cmd_line = ["wopmars", "-w", wopfile_path, "-D", db_url, "-d", self.__testdir_path, "-t", "SampleInformation"]
    #     with self.assertRaises(SystemExit) as se:
    #         WopMars().run(cmd_line)
    #     self.assertEqual(se.exception.code, 1)
    #     #
    #     shutil.rmtree(test_outdir)
    #
    # def test_03sort_reads(self):
    #     """
    #     This test takes as input a db generated by sample_information and the fasta files.
    #     This test write the Variant table
    #     At the end, this test will "delete from Variant"
    #
    #     :return:
    #     """
    #     # input
    #     sample_info_tsv = os.output.join(PathManager.get_module_test_path(), "input", "02sample_information", "fastainfo.tsv")
    #     test_indir = os.output.join(PathManager.get_module_test_path(), "input", "03sort_reads")
    #     fasta_dir = os.output.join(self.__testdir_path, "input", "02sample_information", "fasta")
    #     # output
    #     test_outdir = os.output.join(self.__testdir_path, "output", "03sort_reads")
    #     PathManager.mkdir_p(test_outdir)
    #     # db
    #     db_path = os.output.join(test_indir, "db.sqlite")
    #     db_url = "sqlite:///" + db_path
    #     #
    #     # Create wopfile
    #     template = Template(self.__wopfile_test_str)
    #     wopfile_str = template.render(SAMPLE_INFORMATION_TSV=sample_info_tsv, OUTDIR=test_outdir, FASTA_DIR=fasta_dir)
    #     wopfile_path = os.output.join(test_outdir, "Wopfile_merge.yml")
    #     with open(wopfile_path, 'w') as fout:
    #         self.__wopfile_test_str = fout.write(wopfile_str)
    #     #
    #     cmd_line = ["wopmars", "-w", wopfile_path, "-D", db_url, "-d", self.__testdir_path, "-t", "VariantReadCount"]
    #     p = subprocess.Popen(cmd_line)
    #     p.wait()
    #     #
    #     # Assert Variant table content
    #     output = ('AATTTGAGCAAGATTAATTGGTACATCATTAAGAATAATCATTCGAATTGAATTAAGAACTCCTGGATCATTTATAGGAAATGATCAAATCTATAATTCGATTGTCACTATCCATGCTTTTATTATAATTTTTTTTATAGTTATACCTATTATAATT',)
    #     con = sqlite3.connect(db_path)
    #     cur = con.cursor()
    #     sql = "SELECT sequence from Variant order by sequence;"
    #     cur.execute(sql)
    #     self.assertTrue(list(cur.fetchone()) == list(output))
    #     cur.close()
    #     #
    #     # # Assert variant_count_detail.tsv
    #     # variant_count_detail_tsv = os.output.join(self.__testdir_path, "output", "03sort_reads", "variant_count_detail.tsv")
    #     # variant_count_detail_bak_tsv = os.output.join(self.__testdir_path, "output_bak", "03sort_reads", "variant_count_detail.tsv")
    #     # assert filecmp.cmp(variant_count_detail_tsv, variant_count_detail_bak_tsv)
    #     #
    #     # reset output variant table
    #     cur = con.cursor()
    #     sql = "delete from Variant;"
    #     cur.execute(sql)
    #     cur.close()
    #     con.commit()
    #     con.close()
    #     #
    #     # remove output file
    #     shutil.rmtree(test_outdir)
    #
    # def test_03sort_reads_utilities_read_count(self):
    #     """
    #     """
    #
    # def test_04filter_store_index_below_lfn1_per_replicate(self):
    #     test_outdir = os.output.join(self.__testdir_path, "output", "04filter")
    #     PathManager.mkdir_p(test_outdir)
    #     variant2sample2replicate2count_df_pkl_path = os.output.join(PathManager.get_module_test_path(), "input", "04filter", "variant2sample2replicate2count_df.pkl")
    #     #
    #     # Input
    #     variant2sample2replicate2count_df = pandas.read_pickle(variant2sample2replicate2count_df_pkl_path)
    #     lfn1_per_replicate_threshold = 0.25
    #     #
    #     # Output
    #     variant2sample2replicate2count = Variant2Sample2Replicate2Count(variant2sample2replicate2count_df)
    #     variant2sample2replicate2count.store_index_below_lfn1_per_replicate(lfn1_per_replicate_threshold)
    #     indices_to_drop = variant2sample2replicate2count.indices_to_drop
    #     indices_to_drop_bak = [76, 88, 117, 119]
    #     self.assertTrue(indices_to_drop == indices_to_drop_bak)
    #     #
    #     shutil.rmtree(test_outdir)
    #
    # def test_04filter_store_index_below_lfn2_per_variant(self):
    #     test_outdir = os.output.join(self.__testdir_path, "output", "04filter")
    #     PathManager.mkdir_p(test_outdir)
    #     variant2sample2replicate2count_df_pkl_path = os.output.join(PathManager.get_module_test_path(), "input", "04filter", "variant2sample2replicate2count_df.pkl")
    #     #
    #     # Input
    #     variant2sample2replicate2count_df = pandas.read_pickle(variant2sample2replicate2count_df_pkl_path)
    #     lfn2_per_variant_threshold = 0.025
    #     #
    #     # Output
    #     variant2sample2replicate2count = Variant2Sample2Replicate2Count(variant2sample2replicate2count_df)
    #     variant2sample2replicate2count.store_index_below_lfn2_per_variant(lfn2_per_variant_threshold)
    #     indices_to_drop = variant2sample2replicate2count.indices_to_drop
    #     #
    #     indices_to_drop_bak = [36, 46, 53]
    #     self.assertTrue(indices_to_drop == indices_to_drop_bak)
    #     #
    #     shutil.rmtree(test_outdir)
    #
    # def test_04filter_store_index_below_lfn2_per_replicate_series(self):
    #     test_outdir = os.output.join(self.__testdir_path, "output", "04filter")
    #     PathManager.mkdir_p(test_outdir)
    #     variant2sample2replicate2count_df_pkl_path = os.output.join(PathManager.get_module_test_path(), "input", "04filter", "variant2sample2replicate2count_df.pkl")
    #     #
    #     # Input
    #     variant2sample2replicate2count_df = pandas.read_pickle(variant2sample2replicate2count_df_pkl_path)
    #     lfn2_per_variant_threshold = 0.0075
    #     #
    #     # Output
    #     variant2sample2replicate2count = Variant2Sample2Replicate2Count(variant2sample2replicate2count_df)
    #     variant2sample2replicate2count.store_index_below_lfn2_per_replicate_series(lfn2_per_variant_threshold)
    #     indices_to_drop = variant2sample2replicate2count.indices_to_drop
    #     indices_to_drop_bak = [27, 36, 46, 53, 88, 92, 104, 122, 209]
    #     self.assertTrue(indices_to_drop == indices_to_drop_bak)
    #     #
    #     shutil.rmtree(test_outdir)
    #
    # def test05_taxassign_vsearch(self):
    #     test_outdir = os.output.join(self.__testdir_path, "output", "05taxassign")
