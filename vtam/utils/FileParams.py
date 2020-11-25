import argparse
import os
import sys
import yaml

from vtam.utils.VTAMexception import VTAMexception
from vtam.utils.Logger import Logger
from vtam.utils import constants


class FileParams:

    def __init__(self, params_path=None):

        self.params_path = params_path

        self.params_file_dic = {}
        if not (params_path is None):
            with open(params_path) as fin:
                # The FullLoader parameter handles the conversion from YAML
                # scalar values to Python the dictionary format
                self.params_file_dic = yaml.load(fin, Loader=yaml.SafeLoader) or {}

        self.params_default_dic = constants.get_params_default_dic()


    def is_valid(self):

        """Check if user parameter set is contained in the default parameter set"""
        for k in self.params_file_dic:
            if not (k in self.params_default_dic):
                Logger.instance().error(VTAMexception('Non-valid parameter "{}" in the file "{}"'.format(k, self.params_path)))
                sys.exit(1)
        return True

    def argparse_checker_params_file(self):

        """Returns params_path only if the parameter set is valid"""

        if not os.path.isfile(self.params_path):
            raise argparse.ArgumentTypeError(
                "The file '{}' does not exist. Please fix it.".format(self.params_path))

        if self.is_valid():
            return self.params_path  # return the tsv_path
        else:
            raise argparse.ArgumentTypeError(
                "The format of file '{}' is wrong. Please look at the example in the VTAM "
                "documentation.".format(self.params_path))

    def get_params_dic(self):
        """Return params dic with default or updated values"""

        params_dic = self.params_default_dic.copy()

        for k in self.params_file_dic:
            params_dic[k] = self.params_file_dic[k]

        return params_dic
