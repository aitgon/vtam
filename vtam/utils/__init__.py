import shlex
import subprocess
import sys

from vtam.utils.PathManager import PathManager


def pip_install_vtam_for_tests():

    """This function is used in the tests when the vtam command is run"""

    cmd = '{} -m pip install . -q --upgrade'.format(sys.executable)
    if sys.platform.startswith("win"):
        args = cmd
    else:
        args = shlex.split(cmd)
    subprocess.run(args=args, check=True, cwd=PathManager.instance().get_project_path())
