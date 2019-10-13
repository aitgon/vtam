"""
Example of module documentation which can be
multiple-lined
"""

import errno
import os
import tempfile

from vtam.utils.Singleton import Singleton
from vtam.utils.VTAMexception import VTAMexception


class PathManager(Singleton):

    def __init__(self):
        self.tempdir = None

    def get_tempdir(self):
        """
        Find the Src directory of the project

        :return: the path leading to the src file of the project
        """
        if self.tempdir is None:
            self.tempdir = tempfile.mkdtemp()
        return self.tempdir


    @staticmethod
    def get_wopfile_test_path():
        """
        Find the Src directory of the project

        :return: the path leading to the src file of the project
        """

        wopfile_test_path = os.path.join(os.path.dirname(__file__), "../../test/input/Wopfile_merge.yml")
        return wopfile_test_path

    @staticmethod
    def get_package_path():
        """
        Find the path of the package

        :return: the path leading to the package path
        """

        test_dir_path = os.path.join(os.path.dirname(__file__), "../..")
        return test_dir_path

    @staticmethod
    def get_module_test_path():
        """
        Find the test path of the project

        :return: the path leading to the test path of the project
        """

        test_dir_path = os.path.join(os.path.dirname(__file__), "../../test")
        return test_dir_path


    @staticmethod
    def mkdir_p(path):
        """ Does not fail if directory already exists"""
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


    @staticmethod
    def check_file_exists_and_is_nonempty(path, error_message=None, abspath=False):
        """Checks if file exists and is not empty

        :param error_message: Optional message to help debug the problem
        :param abspath: If True, returns abspath
        :return: void
        """
        try:
            assert os.stat(path).st_size > 0
        except AssertionError as err:
            raise VTAMexception("{}: {}".format(err, error_message))
        if abspath:
            return os.path.abspath(path)
        return path


    @staticmethod
    def check_dir_exists_and_is_nonempty(path, error_message=None, abspath=False):
        """Checks if directory exists and is not empty

        :param error_message: Optional message to help debug the problem
        :param abspath: If True, returns abspath
        :return: void
        """
        try:
            assert len(os.listdir(path)) > 0
            # assert True
        except AssertionError as err:
            raise VTAMexception("{}: {}".format(err, error_message))
        if abspath:
            return os.path.abspath(path)
        return path
