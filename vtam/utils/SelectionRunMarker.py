import argparse

import pandas

from vtam.utils.ArgParser import ArgParserChecker


class SelectionRunMarker(object):

    def __init__(self, path):
        self.path = path
        self.check_argument()

    @staticmethod
    def help():

        return """Input TSV file with two columns and headers 'Run' and 'Marker'.
                                        Example:
                                        Run	Marker
                                        prerun	MFZR
                                        prerun	ZFZR"""

    def check_argument(self):

        """Checks if fastqinfo exists, is not empty

        :param path: Valid non-empty TSV fastqinfo path
        :return: void

        """

        path = ArgParserChecker.check_file_exists_and_is_nonempty(self.path)
        df = pandas.read_csv(path, sep='\t', header=0)
        header = {'Run', 'Marker'}
        if set(df.columns.tolist()) > header:
            raise argparse.ArgumentTypeError("The header of the file {} does not contain these fields: {}!".format(path, header))
        else:
            return path

