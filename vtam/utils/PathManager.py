"""
Example of module documentation which can be
multiple-lined
"""

import os
import pathlib
import tempfile

from vtam.utils.Singleton import Singleton


class PathManager(Singleton):

    def __init__(self):
        self.tempdir = None

    def get_configdir(self):
        """
        Find the ~/.config/vtam dir

        :return: ~/.config/vtam
        """

        vtam_config_dir = os.path.join(str(pathlib.Path.home()), '.config', 'vtam')
        pathlib.Path(vtam_config_dir).mkdir(parents=True, exist_ok=True)

        return vtam_config_dir

    def get_tempdir(self):
        """
        Find the Src directory of the project

        :return: the output leading to the src file of the project
        """
        if self.tempdir is None:
            self.tempdir = tempfile.mkdtemp()
        pathlib.Path(self.tempdir).mkdir(parents=True, exist_ok=True)
        return self.tempdir

    @staticmethod
    def get_package_path():
        """
        Find the output of the precomputed

        :return: the output leading to the precomputed output
        """

        package_path = os.path.join(os.path.dirname(__file__), "../..")
        return package_path

    @staticmethod
    def get_test_path():
        """
        Find the tests output of the project

        :return: the output leading to the tests output of the project
        """

        test_dir_path = os.path.join(os.path.dirname(__file__), "../tests")
        return test_dir_path
