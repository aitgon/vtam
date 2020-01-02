"""
Example of module documentation which can be
multiple-lined
"""

import os
import tempfile

from vtam.utils.Singleton import Singleton


class PathManager(Singleton):

    def __init__(self):
        self.tempdir = None

    def get_tempdir(self):
        """
        Find the Src directory of the project

        :return: the output leading to the src file of the project
        """
        if self.tempdir is None:
            self.tempdir = tempfile.mkdtemp()
        return self.tempdir

    @staticmethod
    def get_package_path():
        """
        Find the output of the precomputed

        :return: the output leading to the precomputed output
        """

        test_dir_path = os.path.join(os.path.dirname(__file__), "../..")
        return test_dir_path

    @staticmethod
    def get_module_test_path():
        """
        Find the test output of the project

        :return: the output leading to the test output of the project
        """

        test_dir_path = os.path.join(os.path.dirname(__file__), "../../test")
        return test_dir_path
