import math

import pandas
import sys

import sqlalchemy
from sqlalchemy import select

from vtam import Logger, VTAMexception
from vtam.utils.SampleInformationUtils import FastaInformationTSV

from vtam.models.Run import Run as RunModel
from vtam.models.Marker import Marker as MarkerModel
from vtam.models.Biosample import Biosample as BiosampleModel
from vtam.models.Variant import Variant
from vtam.models.VariantReadCount import VariantReadCount


class KnownOccurrences(object):

    def __init__(self, known_occurrences_df, readinfo_tsv, engine):
        """
        A class to manipulate the known variant file for the optimize wrappers

        :param known_occurrences_tsv: TSV file with known variants
        """
        # self.known_occurrences_tsv = known_occurrences_tsv
        self.known_occurrences_df = known_occurrences_df
        self.engine = engine
        self.readinfo_tsv = readinfo_tsv

        ################################################################################################################
        #
        # Create a DF and rename columns
        #
        ################################################################################################################

        # self.known_occurrences_df = pandas.read_csv(self.known_occurrences_tsv, sep="\t", header=0, \
        #                                       names=['Marker', 'Run', 'Biosample', 'BiosampleType',
        #                                              'VariantId', 'Action', 'Sequence'], index_col=False,
        #                                         usecols=list(range(7)))
        # Change column names to lower
        self.known_occurrences_df.columns = map(str.lower, self.known_occurrences_df.columns)
        # Rename columns
        self.known_occurrences_df.rename({'marker': 'marker_name',
                                          'run': 'run_name',
                                          'biosample': 'biosample_name',
                                          'biosampletype': 'biosample_type',
                                          'variantid': 'variant_id',
                                          'sequence': 'variant_sequence'}, inplace=True, axis=1)
        # Sequence to upper case
        self.known_occurrences_df['variant_sequence'] = self.known_occurrences_df.variant_sequence.str.upper()

        # columns: run_id, marker_id, biosample_id, replicate, variant_id, biosample_type, action, variant_sequence
        self.known_occurrences_ids_df = pandas.DataFrame.from_records(self.get_ids_of_run_marker_biosample_replicate())

        ################################################################################################################
        #
        # Quality control: Verify if user variant IDs and sequences are consistent with information in the database
        #
        ################################################################################################################

        self.__are_known_variants_coherent_with_db()

        ################################################################################################################
        #
        # Quality control: Verify if run_id, marker_id and biosample_id in known variants are coherent with current analysis (fasta_info)
        #
        ################################################################################################################

        self.__are_known_variants_coherent_with_fasta_info_file()


    def get_ids_of_run_marker_biosample_replicate(self):
        """Returns a list of dictionnaries with run_id, marker_id, biosample_id and replicate entries (See return)

        :return: list of dictionnaries: [{'run_id': 1, 'marker_id': 1, 'biosample_id': 1, 'replicate': 1}, {'run_id': 1, ...
        """
        instance_list = []
        for row in self.known_occurrences_df.itertuples():
            marker_name = row.marker_name
            run_name = row.run_name
            biosample_name = row.biosample_name
            with self.engine.connect() as conn:
                #
                # get run_id ###########
                #
                stmt_select_run_id = sqlalchemy.select([RunModel.__table__.c.id]).where(RunModel.__table__.c.name == run_name)
                run_id_row = conn.execute(stmt_select_run_id).first()
                if run_id_row is None:
                    Logger.instance().error(VTAMexception("Run {} not found in the DB. The program will exit. "
                                                          "Please verify the known_occurrences TSV file.".format(run_name)))
                    sys.exit(1)
                else:
                    run_id = run_id_row[0]
                #
                # get marker_id ###########
                #
                stmt_select_marker_id = sqlalchemy.select([MarkerModel.__table__.c.id]).where(MarkerModel.__table__.c.name == marker_name)
                marker_id_row = conn.execute(stmt_select_marker_id).first()
                if marker_id_row is None:
                    Logger.instance().error(VTAMexception("Marker {} not found in the DB. The program will exit. "
                                                          "Please verify the known_occurrences TSV file.".format(marker_name)))
                    sys.exit(1)
                else:
                    marker_id = marker_id_row[0]
                #
                # get biosample_id ###########
                #
                stmt_select_biosample_id = sqlalchemy.select([BiosampleModel.__table__.c.id]).where(BiosampleModel.__table__.c.name == biosample_name)
                biosample_id_first = conn.execute(stmt_select_biosample_id).first()
                if biosample_id_first is None:
                    Logger.instance().error(VTAMexception("Biosample {} not found in the DB. The program will exit. "
                                                          "Please verify the known_occurrences TSV file.".format(biosample_name)))
                    sys.exit(1)
                else:
                    biosample_id = biosample_id_first[0]
                # get replicate ###########
                # add this sample_instance ###########
                instance_list.append({'run_id': run_id, 'marker_id': marker_id, 'biosample_id': biosample_id, 'variant_id': row.variant_id,
                                        'biosample_type': row.biosample_type, 'action': row.action, 'variant_sequence': row.variant_sequence})
        return instance_list


    def get_known_occurrences_df(self):
        """

        :return: known_occurrences in variant_read_count_input_df format
        """
        return self.known_occurrences_df

    def __are_known_variants_coherent_with_db(self):
        """Raises an error and exists if some of these conditions
            - the sequence is not None (negative biosample) and does not exist in the DB
            - the variant ID/sequence are not coherent with the database

        :return void
        """

        variant_control_df = self.known_occurrences_df[['variant_id', 'variant_sequence']].drop_duplicates()
        variant_control_df = variant_control_df.loc[~variant_control_df.variant_sequence.isnull()]

        with self.engine.connect() as conn:

            # for row in variant_control_df.itertuples():
            for row in self.known_occurrences_ids_df.itertuples():

                # import pdb; pdb.set_trace()

                ###########################################################################
                #
                # Checks if variant sequence exists in the DB
                #
                ###########################################################################

                user_variant_id = row.variant_id
                user_variant_sequence = row.variant_sequence
                user_run_id = row.run_id
                user_marker_id = row.marker_id
                user_biosample_id = row.biosample_id

                # Implicit 'delete' variants in 'negative' biosamples do not have sequence
                if not (user_variant_sequence is None):  # variants have sequence

                    stmt_select = select([VariantReadCount.__table__.c.variant_id, Variant.__table__.c.sequence])\
                        .where(Variant.__table__.c.sequence == user_variant_sequence)\
                        .where(VariantReadCount.__table__.c.variant_id == Variant.__table__.c.id)\
                        .where(VariantReadCount.__table__.c.run_id == user_run_id)\
                        .where(VariantReadCount.__table__.c.marker_id == user_marker_id)\
                        .where(VariantReadCount.__table__.c.biosample_id == user_biosample_id)

                    known_variant_in_db = conn.execute(stmt_select).first()

                    if known_variant_in_db is None:  # user sequence not found in db

                        msg_error = 'Error: This variant sequence was not found in the given run-marker-biosample combination in the DB. ' \
                                    'Please verify this variant sequence: {}'\
                            .format(user_variant_sequence)
                        Logger.instance().error(VTAMexception(msg_error))
                        sys.exit(1)

                    db_variant_id, db_variant_sequence = known_variant_in_db

                    if math.isnan(user_variant_id):  # User has not given variant id, then keep this variant id

                        self.known_occurrences_df.loc[
                            self.known_occurrences_df.variant_sequence == user_variant_sequence, 'variant_id'] = db_variant_id

                        self.known_occurrences_ids_df.loc[
                            self.known_occurrences_ids_df.variant_sequence == user_variant_sequence, 'variant_id'] = db_variant_id

                    ###########################################################################
                    #
                    # Checks if variant id/sequence is coherent in the DB
                    #
                    ###########################################################################

                    else:  # user has given variant id, check if coherent

                        if not user_variant_id == db_variant_id: # user variant id does not correspond to sequence in DB

                            msg_error = 'Error: This combination of variant id / sequence was not found in the DB. ' \
                                    'Please verify the variant ID for the given sequence: {}'\
                                    .format(user_variant_id, user_variant_sequence)
                            Logger.instance().error(VTAMexception(msg_error))
                            sys.exit(1)

    def __are_known_variants_coherent_with_fasta_info_file(self):

        fasta_info_tsv = FastaInformationTSV(engine=self.engine, fasta_info_tsv=self.readinfo_tsv)
        sample_information_df = fasta_info_tsv.sample_information_df

        ################################################################################################################
        #
        # Takes fasta_info and intersects with known variants to create a variant known variant_read_count_input_df of this experiment
        #
        ################################################################################################################

        for row in self.known_occurrences_ids_df.itertuples():
            run_id = row.run_id
            marker_id = row.marker_id
            biosample_id = row.biosample_id
            try:
                assert (sample_information_df.loc[
                    (sample_information_df['biosample_id'] == biosample_id) & (sample_information_df['marker_id'] == marker_id) & (
                                sample_information_df['run_id'] == run_id)]).shape[0] > 0
            except AssertionError:
                Logger.instance().error(VTAMexception("Error: Verify in the --known_occurrences file that run_id, marker_id and biosample_id"
                                                      "are defined in the --fasta_info file"))
                sys.exit(1)

    def get_keep_run_marker_biosample_variant_df(self, variant_tolerate=False):
        """Returns the 'keep' and 'tolerates' variants together with run_id, marker_id, biosample_id, variant_id


        :param: variant_tolerate: Boolean: Default False. include "variant_tolerate" variants or not?
        :return: pandas variant_read_count_input_df with columns: run_id, marker_id, biosample_id, variant_id
        """
        # Get portion of known_occurrences_tsv with either keep or keep+variant_tolerate
        if variant_tolerate:  # get also variant_tolerate variant
            run_marker_biosample_variant_keep_df = self.known_occurrences_ids_df.loc[
                ((self.known_occurrences_df.action.isin(['keep', 'variant_tolerate']))).values]
        else:  # do only get keep variants
            run_marker_biosample_variant_keep_df = self.known_occurrences_ids_df.loc[(self.known_occurrences_df.action == 'keep').values]
        # run_marker_biosample_variant_keep_df = run_marker_biosample_variant_keep_df.merge(self.known_occurrences_ids_df,
        #                                         on=['run_id', 'marker_id', 'biosample_id', 'variant_id'])
        # Select run_id, marker_id, biosample_id and variant_id
        run_marker_biosample_variant_keep_df = run_marker_biosample_variant_keep_df[
            ['run_id', 'marker_id', 'biosample_id', 'variant_id']].drop_duplicates(inplace=False)
        # Change variant_id type to int
        # run_marker_biosample_variant_keep_df.variant_id = run_marker_biosample_variant_keep_df.variant_id.astype('int')
        return run_marker_biosample_variant_keep_df


    def get_delete_run_marker_biosample_variant_df(self, variant_read_count_df):
        """Returns the 'delete' variants together with run_id, marker_id, biosample_id, variant_id

        :return: pandas variant_read_count_input_df with columns: run_id, marker_id, biosample_id, variant_id
        """

        ##########################################################
        #
        # Get delete variants, that are not keep and/or tolerate in mock or real
        #
        ##########################################################
        # Get mock biosamples
        run_marker_biosample_mock_df = self.known_occurrences_ids_df.loc[
            (self.known_occurrences_df.biosample_type == 'mock').values, ['run_id', 'marker_id', 'biosample_id']]
        # Get variant_read_count_mock
        variant_read_count_mock = run_marker_biosample_mock_df.merge(variant_read_count_df,
                                                                     on=['run_id', 'marker_id', 'biosample_id'])
        # Get keep and tolerates variants to remove from mock
        run_marker_biosample_variant_keep_df = self.get_keep_run_marker_biosample_variant_df(variant_tolerate=True)
        # Compute variants that are 1) mock and keep/tolerate (both). 2) Mock but not keep/tolerate (left_only)
        variant_delete_mock_df = variant_read_count_mock.merge(run_marker_biosample_variant_keep_df,
                                   on=['run_id', 'marker_id', 'biosample_id', 'variant_id'], how='left', indicator=True)
        # Throw replicate and read count
        variant_delete_mock_df = variant_delete_mock_df.loc[variant_delete_mock_df._merge=='left_only',
                                    ['run_id', 'marker_id', 'biosample_id', 'variant_id']].drop_duplicates(inplace=False)

        ##########################################################
        #
        # Get delete variants, that are in negative samples
        #
        ##########################################################
        # Get negative biosamples
        variant_delete_negative_df = self.known_occurrences_ids_df.loc[(self.known_occurrences_df.biosample_type == 'negative').values,
                                                          ['run_id', 'marker_id', 'biosample_id']]
        # Inner merge of variants and negative biosamples
        variant_delete_negative_df = variant_delete_negative_df.merge(variant_read_count_df, on=['run_id', 'marker_id',
                                                                                               'biosample_id'])
        # Throw replicate and read count
        variant_delete_negative_df = variant_delete_negative_df[['run_id', 'marker_id', 'biosample_id', 'variant_id']] \
            .drop_duplicates(inplace=False)

        ##########################################################
        #
        # Get delete variants, that are marked so in any (real) samples
        #
        ##########################################################
        # Get run_id, marker_id, biosample_id, variant_id that are explicitely marked as delete
        variant_delete_real_df = self.known_occurrences_ids_df.loc[(self.known_occurrences_df.action == 'delete').values]
        # Remove delete variant that are not explicite, ie in negative biosamples
        variant_delete_real_df = variant_delete_real_df[~variant_delete_real_df.variant_id.isnull()]
        #  Throw replicate and read count
        variant_delete_real_df = variant_delete_real_df[
            ['run_id', 'marker_id', 'biosample_id', 'variant_id']].drop_duplicates(inplace=False)

        ##########################################################
        #
        # Merge (Vertically) the three classes of delete variants
        #
        ##########################################################
        variant_delete_df = pandas.concat([variant_delete_mock_df, variant_delete_negative_df, variant_delete_real_df])
        variant_delete_df = variant_delete_df.drop_duplicates(inplace=False)
        variant_delete_df.variant_id = variant_delete_df.variant_id.astype(int)
        variant_delete_df = variant_delete_df.reset_index(drop=True)
        variant_delete_df.drop_duplicates(inplace=True)
        return variant_delete_df, variant_delete_mock_df, variant_delete_negative_df, variant_delete_real_df
