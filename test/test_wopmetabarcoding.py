import filecmp
import os
import sqlite3
from unittest import TestCase
import subprocess

from wopmetabarcoding.utils.PathFinder import PathFinder


class TestWopMetabarcoding(TestCase):
    def setUp(self):
        self.__db_path = os.path.join(PathFinder.get_module_test_path(), "output/db.sqlite")
        self.__db_url = "sqlite:///" + self.__db_path

    def test_01run(self):
        sample_info_csv = os.path.join(PathFinder.get_module_test_path(), "input", "1_sample_info.csv")
        cmd_line = ["wopmars", "-w", PathFinder.get_wopfile_path(), "-D", self.__db_url, "-d", PathFinder.get_module_path(), "-p", "-v", "-F", "-t", "ReadCSV"]
        p = subprocess.Popen(cmd_line)
        p.wait()
        #
        backup_table_list = [('wom_execution',), ('wom_modification_table',), ('wom_type',), ('wom_rule',), ('wom_option',), ('wom_table',), ('wom_file',), ('Marker',), ('PrimerPair',), ('Variant',), ('TagPair',), ('Replicate',), ('Replicatemarker',), ('SampleInformation',), ('File',), ('ReadCount',), ('Biosample',)]
        con = sqlite3.connect(self.__db_path)
        cur = con.cursor()
        cur.execute("SELECT * from SampleInformation;")
        print(cur.fetchall())
        self.assertTrue(cur.fetchall() == backup_table_list)
