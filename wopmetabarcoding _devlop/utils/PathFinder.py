"""
Example of module documentation which can be
multiple-lined
"""

import errno
import os


class PathFinder:


    """
    Static class for finding paths
    """
    @staticmethod
    def get_module_path():
        """
        Find the Src directory of the project

        :return: the path leading to the src file of the project
        """

        module_path = os.path.join(os.path.dirname(__file__), "../..")
        return module_path

    @staticmethod
    def get_wopfile_test_path():
        """
        Find the Src directory of the project

        :return: the path leading to the src file of the project
        """

        wopfile_path = os.path.join(os.path.dirname(__file__), "../../test/input/Wopfile.yml")
        return wopfile_path

    @staticmethod
    def get_module_test_path():
        """
        Find the Src directory of the project

        :return: the path leading to the src file of the project
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

