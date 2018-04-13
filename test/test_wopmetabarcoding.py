import filecmp
import os
import shutil
import sqlite3
from unittest import TestCase
import subprocess
from jinja2 import Template

from wopmetabarcoding.utils.PathFinder import PathFinder



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
        test_outdir = os.path.join(self.__testdir_path, "output", "sample_information")
        PathFinder.mkdir(test_outdir)
        sample_info_csv = os.path.join(PathFinder.get_module_test_path(), "input", "sample_information", "1_sample_info.csv")
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
        cmd_line = ["wopmars", "-w", wopfile_path, "-D", db_url, "-d", PathFinder.get_module_path(), "-F", "-t", "SampleInformation"]
        p = subprocess.Popen(cmd_line)
        p.wait()
        #
        output = (1, 'ZFZR', 'AGATATTGGAACWTTATATTTTATTTTTGG', 'WACTAATCAATTWCCAAATCCTCC', 'cgatcgtcatcacg', 'ctcgatgatcacg', 1, 'prerun', '14Mon01', 'repl1')
        print(self.__db_path)
        con = sqlite3.connect(self.__db_path)
        print(con)
        cur = con.cursor()
        cur.execute("SELECT * from SampleInformation;")
        self.assertTrue(list(cur.fetchone()) == list(output))
        #
        shutil.rmtree(test_outdir)

    def test_02_sort_reads(self):
        test_outdir = os.path.join(self.__testdir_path, "output", "sort_reads")
        PathFinder.mkdir(test_outdir)
        sample_info_csv = os.path.join(PathFinder.get_module_test_path(), "input", "sort_reads", "1_sample_info.csv")
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
        marker = 'ZFZR'
        sequence = 'AGCCTGAGCTGGAATAGTAGGTACTTCCCTTAGTATACTTATTCGAGCCGAATTAGGACACCCAGGCTCTCTAATTGGAGACGACCAAATTTATAATGTAATTGTTACTGCTCATGCTTTTGTAATAATTTTTTTTATAGTTATGCCAATTATAATT'
        con = sqlite3.connect(db_path)
        cur = con.cursor()
        sql = "SELECT marker from Variant where sequence='%s';"%(sequence)
        cur.execute(sql)
        self.assertTrue(cur.fetchone()[0] == marker)
        #
        shutil.rmtree(test_outdir)

