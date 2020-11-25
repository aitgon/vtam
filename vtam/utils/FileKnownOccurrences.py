import argparse
import os
import pandas

from vtam.models.Sample import Sample
from vtam.models.Marker import Marker
from vtam.models.Run import Run
from vtam.utils.constants import header_known_occurrences

from vtam.utils.NameIdConverter import NameIdConverter


class FileKnownOccurrences(object):

    def __init__(self, known_occurrences_tsv):
        """
        A class to manipulate the known variant file for the optimize wrappers

        :param tsv_path: TSV file with known variants
        """
        self.known_occurrences_tsv = known_occurrences_tsv

    def read_tsv_into_df(self):
        """Read into stats_df
        Updated: Mai 17, 2020

        Parameters
        ----------

        Returns
        -------
        pandas.DataFrame

        """

        known_occurrences_df = pandas.read_csv(
            self.known_occurrences_tsv, sep="\t", header=0)
        known_occurrences_df.columns = known_occurrences_df.columns.str.lower()
        known_occurrences_df = known_occurrences_df[header_known_occurrences]
        known_occurrences_df.rename({'run': 'run_name', 'marker': 'marker_name',
                                     'sample': 'sample_name', 'variant': 'variant_id',
                                     'sequence': 'variant_sequence'}, axis=1,
                                    inplace=True)
        return known_occurrences_df

    def argparse_checker_known_occurrences(self):

        """
        Used by the argparser to check whether the file is in the correct format
        Updated: June 2, 2020

        """

        if not os.path.isfile(self.known_occurrences_tsv):
            raise argparse.ArgumentTypeError(
                "The file '{}' does not exist. Please fix it.".format(self.known_occurrences_tsv))
        elif not os.stat(self.known_occurrences_tsv).st_size > 0:
            raise argparse.ArgumentTypeError(
                "The file '{}' is empty!".format(self.known_occurrences_tsv))
        known_occurrences_df = pandas.read_csv(
            self.known_occurrences_tsv, sep="\t", header=0)
        known_occurrences_df.columns = known_occurrences_df.columns.str.lower()
        # must contain at least these columns
        if set(known_occurrences_df.columns) >= header_known_occurrences:
            return self.known_occurrences_tsv  # return the tsv_path
        else:
            raise argparse.ArgumentTypeError(
                "The format of file '{}' is wrong. Please look at the example in the VTAM "
                "documentation.".format(self.known_occurrences_tsv))

    def to_identifier_df(self, engine):
        """Returns a list of dictionnaries with run_id, marker_id, sample_id entries (See return)

        :return: pandas.DataFrame: with columns run_id, marker_id, ...
        """
        instance_list = []
        df = self.read_tsv_into_df()

        df.run_name = NameIdConverter(df.run_name.tolist(), engine).to_ids(Run)
        df.marker_name = NameIdConverter(df.marker_name.tolist(), engine).to_ids(Marker)
        df.sample_name = NameIdConverter(df.sample_name.tolist(), engine).to_ids(Sample)
        df['variant_id'] = NameIdConverter(df.variant_sequence.tolist(), engine)\
            .variant_sequence_to_id()
        df.rename({'run_name': 'run_id', 'marker_name': 'marker_id',
                   'sample_name': 'sample_id'}, axis=1, inplace=True)
        return df

    def get_run_marker_sample_variant_df(self, engine, action):
        """

        :param engine: sqlalchemy engine necessary to ID data
        :param action: takes values 'keep' or 'delete'
        :return:
        """

        """Returns the 'keep' and 'tolerates' variants together with run_id, marker_id, sample_id, variant_id

        :param: variant_tolerate: Boolean: Default False. include "variant_tolerate" variants or not?
        :return: pandas variant_read_count_input_df with columns: run_id, marker_id, sample_id, variant_id
        """

        if not (action in ['keep', 'delete']):
            raise argparse.ArgumentTypeError(
                "The argument 'action' takes either 'keep' or 'delete' values")

        # Get portion of tsv_path with either keep or delete
        run_marker_sample_variant_keep_df = self.to_identifier_df(engine).loc[(
            (self.read_tsv_into_df()).action == action).values]
        # Select run_id, marker_id, sample_id and variant_id
        run_marker_sample_variant_keep_df = run_marker_sample_variant_keep_df[[
            'run_id', 'marker_id', 'sample_id', 'variant_id']].drop_duplicates(inplace=False)
        # Change variant_id type to int
        # run_marker_sample_variant_keep_df.variant_id = run_marker_sample_variant_keep_df.variant_id.astype('int')
        return run_marker_sample_variant_keep_df
