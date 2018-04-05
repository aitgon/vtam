import filecmp
import os
import sqlite3
from unittest import TestCase
import subprocess
import jinja2

from wopmetabarcoding.utils.PathFinder import PathFinder


class TestWopMetabarcoding(TestCase):
    def setUp(self):
        self.__db_path = os.path.join(PathFinder.get_module_test_path(), "output/db.sqlite")
        self.__db_url = "sqlite:///" + self.__db_path
        self.__wopfile_template_path = PathFinder.get_wopfile_template_path()
        self.__wopfile_template_str = ""
        with open(self.__wopfile_template_path, 'r') as fin:
            self.__wopfile_template_str = fin.read()

    def test_01_ReadCSV(self):
        sample_info_csv = os.path.join(PathFinder.get_module_test_path(), "input", "1_sample_info.csv")
        #
        from jinja2 import Template
        template = Template(self.__wopfile_template_str)
        # print(template.render(SAMPLE_INFORMATION_CSV=sample_info_csv))
        wopfile_path = os.path.join(PathFinder.get_module_test_path(), "output", "Wopfile.yml")
        with open(wopfile_path, 'w') as fout:
            self.__wopfile_template_str = fout.write(template.render(SAMPLE_INFORMATION_CSV=sample_info_csv))
        #
        cmd_line = ["wopmars", "-w", wopfile_path, "-D", self.__db_url, "-d", PathFinder.get_module_path(), "-F", "-t", "ReadCSV"]
        p = subprocess.Popen(cmd_line)
        p.wait()
        #
        backup_table_list = [(1, 'MFZR', 'TCCACTAATCACAARGATATTGGTAC', 'WACTAATCAATTWCCAAATCCTCC', 'cgatcgtcatcacg', 'cacgatttgtagag', 'data/fastq_merged_fasta/prerun_MFZR_repl1.fasta', 'prerun', '14Mon01', 'repl1'), (2, 'MFZR', 'TCCACTAATCACAARGATATTGGTAC', 'WACTAATCAATTWCCAAATCCTCC', 'tcgatcacgatgt', 'cgcgtctgtagag', 'data/fastq_merged_fasta/prerun_MFZR_repl1.fasta', 'prerun', '14Mon02', 'repl1')]
        con = sqlite3.connect(self.__db_path)
        cur = con.cursor()
        cur.execute("SELECT * from SampleInformation;")
        self.assertTrue(cur.fetchall() == backup_table_list)
