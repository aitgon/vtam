import argparse
import os

import pandas

from vtam.utils.constants import header_cutoff_specific_variant_replicate, \
    header_cutoff_specific_variant


class CutoffSpecificFile(object):

    def __init__(self, cutoff_specific_tsv):
        """
        A class to manipulate the known variant file for the optimize wrappers

        :param tsv_path: TSV file with known variants
        """
        self.cutoff_specific_tsv = cutoff_specific_tsv

    def argparse_checker(self):

        """
        Used by the argparser to check whether the file is in the correct format
        Updated: June 2, 2020

        """

        if not os.path.isfile(self.cutoff_specific_tsv):
            raise argparse.ArgumentTypeError(
                "The file '{}' does not exist. Please fix it.".format(self.cutoff_specific_tsv))
        elif not os.stat(self.cutoff_specific_tsv).st_size > 0:
            raise argparse.ArgumentTypeError(
                "The file '{}' is empty!".format(self.cutoff_specific_tsv))
        cutoff_specific_df = pandas.read_csv(
            self.cutoff_specific_tsv, sep="\t", header=0)
        cutoff_specific_df.columns = cutoff_specific_df.columns.str.lower()
        # must contain at least these columns
        if (set(cutoff_specific_df.columns) >= header_cutoff_specific_variant_replicate) \
                or (set(cutoff_specific_df.columns) >= header_cutoff_specific_variant):
            return self.cutoff_specific_tsv  # return the tsv_path
        else:
            raise argparse.ArgumentTypeError(
                "The format of file '{}' is wrong. Probably the file does not contain this column "
                "header{} or {}. Please look at the example in the VTAM  documentation."
                    .format(self.cutoff_specific_tsv, header_cutoff_specific_variant,
                            header_cutoff_specific_variant_replicate))
