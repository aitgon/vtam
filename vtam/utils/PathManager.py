"""
Example of module documentation which can be
multiple-lined
"""

import os
import pathlib
import tempfile
import vtam

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

    @classmethod
    def get_doc_path(cls):
        """
        Returns the document folder

        :return: path to the document folder
        """

        doc_path = os.path.join(cls.get_package_path(), "../doc")
        return doc_path

    @classmethod
    def get_project_path(cls):
        """
        Returns the path to the project root

        :return: path to the root of the project
        """

        project_path = os.path.join(cls.get_package_path(), "..")
        return project_path

    @staticmethod
    def get_package_path():
        """
        Returns the vtam.__path__[0]

        :return: path to the package
        """

        package_path = vtam.__path__[0]
        return package_path

    @classmethod
    def get_test_path(cls):
        """
        Find the tests output of the project

        :return: the output leading to the tests output of the project
        """

        test_path = os.path.join(cls.get_package_path(), "tests")
        return test_path
