"""
Example of module documentation which can be
multiple-lined
"""

import errno
import os

from wopmetabarcoding.utils.VTAMexception import VTAMexception


class PathFinder:

    @staticmethod
    def get_wopfile_test_path():
        """
        Find the Src directory of the project

        :return: the path leading to the src file of the project
        """

        wopfile_path = os.path.join(os.path.dirname(__file__), "../../test/input/Wopfile_merge.yml")
        return wopfile_path

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
    def check_file_exists_and_is_nonempty(path, error_message=None):
        """Checks if file exists and is not empty

        :param error_message: Optional message to help debug the problem
        :return: void
        """
        try:
            if os.stat(path).st_size > 0:
                return os.path.abspath(path)
        except OSError as err:
            raise VTAMexception("{}: {}".format(err, error_message))


    @staticmethod
    def check_dir_exists_and_is_nonempty(path, error_message=None):
        """Checks if directory exists and is not empty

        :param error_message: Optional message to help debug the problem
        :return: void
        """
        try:
            assert len(os.listdir(path)) > 0
            # assert True
        except NotADirectoryError as err:
            raise VTAMexception("{}: {}".format(err, error_message))
        return os.path.abspath(path)

