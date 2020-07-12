from vtam.utils.PathManager import PathManager
import os
import pathlib
import shlex
import shutil
import sqlite3
import subprocess
import sys
import unittest


class TestCmdVariantReadCount(unittest.TestCase):

    """Will test main commands based on a complete test dataset"""

    @classmethod
    def setUpClass(cls):

        # vtam needs to be in the tsv_path
        cmd = '{} -m pip install . -q --upgrade'.format(sys.executable)
        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=PathManager.get_package_path())

    def setUp(self):

        self.package_path = os.path.join(PathManager.get_package_path())
        self.test_path = os.path.join(PathManager.get_test_path())
        self.outdir_path = os.path.join(self.test_path, 'outdir')
        shutil.rmtree(self.outdir_path, ignore_errors=True)
        pathlib.Path(self.outdir_path).mkdir(parents=True, exist_ok=True)

    def test_filter_singleton(self):

        args = {}
        args['readinfo'] = os.path.join(os.path.dirname(__file__), "readinfo_1sample_singletons.tsv")
        args['sorteddir'] = os.path.dirname(__file__)

        cmd = "vtam filter --db db.sqlite --readinfo {readinfo} --readdir {sorteddir} " \
              "--asvtable asvtable_default.tsv --until VariantReadCount".format(**args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        result = subprocess.run(args=args, cwd=self.outdir_path)

        self.assertEqual(result.returncode, 0)

    def test_filter_fasta_over_two_lines(self):

        args = {}
        args['readinfo'] = os.path.join(os.path.dirname(__file__), "readinfo_1sample.tsv")
        args['sorteddir'] = os.path.dirname(__file__)

        cmd = "vtam filter --db db.sqlite --readinfo {readinfo} --readdir {sorteddir} " \
              "--asvtable asvtable_default.tsv --until VariantReadCount".format(**args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        result = subprocess.run(args=args, cwd=self.outdir_path)

        self.assertEqual(result.returncode, 0)

        db_path = os.path.join(self.outdir_path, "db.sqlite")
        con = sqlite3.connect(db_path)
        cur = con.cursor()
        cur_result = cur.execute('SELECT COUNT(*) from Variant').fetchone()
        self.assertEqual(cur_result[0], 2)
        con.close()

    def tearDown(self):

        shutil.rmtree(self.outdir_path, ignore_errors=True)
