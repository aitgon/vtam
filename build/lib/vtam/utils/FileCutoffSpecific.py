import argparse
import os
import sys
import pandas

from vtam.utils.Logger import Logger
from vtam.utils.VTAMexception import VTAMexception
from vtam.models.Marker import Marker
from vtam.models.Run import Run
from vtam.utils.NameIdConverter import NameIdConverter
from vtam.utils.constants import header_cutoff_specific_variant_replicate, \
    header_cutoff_specific_variant


class FileCutoffSpecific(object):

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

    def is_compatible_lfn_variant_replicate(self):

        """
        Verifies if this cutoff specific files contains header_cutoff_specific_variant_replicate columns
        Updated: June 18, 2020

        """


        cutoff_specific_df = pandas.read_csv(
            self.cutoff_specific_tsv, sep="\t", header=0)
        cutoff_specific_df.columns = cutoff_specific_df.columns.str.lower()
        # must contain at least these columns
        if (set(cutoff_specific_df.columns) >= header_cutoff_specific_variant_replicate):
            return True
        else:
            return False

    def read_tsv_into_df(self, is_lfn_variant_replicate):
        """Read into stats_df
        Updated: June 3, 2020

        Parameters
        ----------
        is_lfn_variant_replicate: Boolean that tells if algorithm is lfn_variant_replicate or not

        Returns
        -------
        pandas.DataFrame

        """

        df = pandas.read_csv(self.cutoff_specific_tsv, sep="\t", header=0)
        df.columns = df.columns.str.lower()
        df.rename({'lfn_variant_cutoff': 'cutoff', 'lfn_variant_replicate_cutoff': 'cutoff'}, inplace=True, axis=1)

        if is_lfn_variant_replicate and set(df.columns.tolist()) >= {'run', 'marker', 'variant', 'replicate', 'cutoff', 'sequence'}:
            df = df[['run', 'marker', 'variant', 'replicate', 'cutoff', 'sequence']]
        elif not is_lfn_variant_replicate and set(df.columns.tolist()) >= {'run', 'marker', 'variant', 'cutoff', 'sequence'}:
            df = df[['run', 'marker', 'variant', 'cutoff', 'sequence']]
        else:
            Logger.instance().critical(VTAMexception("The format of file '{}' is wrong. Columns 'lfn_variant_cutoff' or 'lfn_variant_replicate_cutoff' are required."
                    .format(self.cutoff_specific_tsv)))
            sys.exit(1)

        df.rename({'run': 'run_name', 'marker': 'marker_name', 'variant': 'variant_id', 'sequence': 'variant_sequence'}, axis=1, inplace=True)

        return df

    def to_identifier_df(self, engine, is_lfn_variant_replicate):
        """Returns a list of dictionnaries with run_id, marker_id, sample_id entries (See return)

        :return: pandas.DataFrame: with columns run_id, marker_id, ...
        """

        df = self.read_tsv_into_df(is_lfn_variant_replicate)

        df.run_name = NameIdConverter(df.run_name.tolist(), engine).to_ids(Run)
        df.marker_name = NameIdConverter(df.marker_name.tolist(), engine).to_ids(Marker)

        variant_id_user_lst = df.variant_id.tolist()
        df['variant_id'] = NameIdConverter(df.variant_sequence.tolist(), engine).variant_sequence_to_id()
        if not df.variant_id.tolist() == variant_id_user_lst:
            Logger.instance().warning(VTAMexception("Some variant IDs and sequences do not agree in the --cutoff_specific file and in the database."))

        df.rename({'run_name': 'run_id', 'marker_name': 'marker_id'}, axis=1, inplace=True)

        return df
