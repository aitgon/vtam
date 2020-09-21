import os
import pathlib
import shlex
import shutil
import subprocess
import sys
import unittest

from vtam.utils.PathManager import PathManager


class TestCommandExample(unittest.TestCase):

    """Will test main commands based on a complete test dataset"""

    def setUp(self):

        self.test_path = PathManager.get_test_path()
        self.outdir_path = os.path.join(self.test_path, 'outdir')
        pathlib.Path(self.outdir_path).mkdir(exist_ok=True, parents=True)

    def test_command_example(self):

        cmd = "vtam example"

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)

        result = subprocess.run(args=args, check=True, cwd=self.outdir_path)
        self.assertEqual(result.returncode, 0)

    def tearDown(self):
        shutil.rmtree(self.outdir_path, ignore_errors=True)
