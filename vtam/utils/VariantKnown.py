import pandas
import sys

import sqlalchemy
from sqlalchemy import select

from vtam import Logger, VTAMexception
from vtam.utils.SampleInformationId import FastaInformation2


class VariantKnown(object):

    def __init__(self, variant_known_tsv, fasta_info_tsv, engine, variant_model, run_model, marker_model, biosample_model):
        """
        A class to manipulate the known variant file for the optimize wrappers

        :param variant_known_tsv: TSV file with known variants
        """
        self.variant_known_tsv = variant_known_tsv
        self.engine = engine
        self.variant_model = variant_model
        self.run_model = run_model
        self.marker_model = marker_model
        self.biosample_model = biosample_model
        self.fasta_info_tsv = fasta_info_tsv

        ################################################################################################################
        #
        # Create a DF
        #
        ################################################################################################################
        self.variant_known_df = pandas.read_csv(self.variant_known_tsv, sep="\t", header=0, \
                                              names=['marker_name', 'run_name', 'biosample_name', 'biosample_type',
                                                     'variant_id', 'action', 'variant_sequence', 'note'], index_col=False)
        # columns: run_id, marker_id, biosample_id, replicate, variant_id, biosample_type, action, variant_sequence
        self.variant_known_ids_df = pandas.DataFrame.from_records(self.get_ids_of_run_marker_biosample_replicate())

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
        for row in self.variant_known_df.itertuples():
            marker_name = row.marker_name
            run_name = row.run_name
            biosample_name = row.biosample_name
            with self.engine.connect() as conn:
                # get run_id ###########
                stmt_select_run_id = sqlalchemy.select([self.run_model.__table__.c.id]).where(self.run_model.__table__.c.name==run_name)
                run_id = conn.execute(stmt_select_run_id).first()[0]
                # get marker_id ###########
                stmt_select_marker_id = sqlalchemy.select([self.marker_model.__table__.c.id]).where(self.marker_model.__table__.c.name==marker_name)
                marker_id = conn.execute(stmt_select_marker_id).first()[0]
                # get biosample_id ###########
                stmt_select_biosample_id = sqlalchemy.select([self.biosample_model.__table__.c.id]).where(self.biosample_model.__table__.c.name==biosample_name)
                biosample_id = conn.execute(stmt_select_biosample_id).first()[0]
                # get replicate ###########
                # add this sample_instance ###########
                instance_list.append({'run_id': run_id, 'marker_id': marker_id, 'biosample_id': biosample_id, 'variant_id': row.variant_id,
                                        'biosample_type': row.biosample_type, 'action': row.action, 'variant_sequence': row.variant_sequence})
        return instance_list


    def get_variant_known_df(self):
        """

        :return: variant_known in df format
        """
        return self.variant_known_df

    def __are_known_variants_coherent_with_db(self):
        """Raises an error and exists if some of the variant ID/sequence are not coherent with the database

        :return void
        """

        variant_control_df = self.variant_known_df[['variant_id', 'variant_sequence']].drop_duplicates()
        variant_control_df = variant_control_df.loc[~variant_control_df.variant_id.isnull()]
        with self.engine.connect() as conn:
            for row in variant_control_df.itertuples():
                variant_id = row.variant_id
                variant_sequence = row.variant_sequence
                stmt_select = select([self.variant_model.__table__.c.id, self.variant_model.__table__.c.sequence])\
                    .where(self.variant_model.__table__.c.id == variant_id)\
                    .where(self.variant_model.__table__.c.sequence == variant_sequence)
                try:
                    assert not conn.execute(stmt_select).first() is None
                except AssertionError:
                    Logger.instance().error(VTAMexception("Error: The variant_id {} and variant_sequence {} from the "
                                                          "known variant file do not match in the DB".format(str(variant_id), variant_sequence)))
                    sys.exit(1)



    def __are_known_variants_coherent_with_fasta_info_file(self):

        fasta_info_obj = FastaInformation2(engine=self.engine, fasta_info_tsv=self.fasta_info_tsv, run_model=self.run_model,
                                      marker_model=self.marker_model, biosample_model=self.biosample_model)
        sample_information_id_df = fasta_info_obj.sample_information_id_df

        ################################################################################################################
        #
        # Takes fasta_info and intersects with known variants to create a variant known df of this experiment
        #
        ################################################################################################################

        for row in self.variant_known_ids_df.itertuples():
            run_id = row.run_id
            marker_id = row.marker_id
            biosample_id = row.biosample_id
            try:
                assert (sample_information_id_df.loc[
                    (sample_information_id_df['biosample_id'] == biosample_id) & (sample_information_id_df['marker_id'] == marker_id) & (
                                sample_information_id_df['run_id'] == run_id)]).shape[0] > 0
            except AssertionError:
                Logger.instance().error(VTAMexception("Error: Verify in the --variant_known file that run_id, marker_id and biosample_id"
                                                      "are defined in the --fasta_info file"))
                sys.exit(1)

    def get_keep_run_marker_biosample_variant_df(self, variant_tolerate=False):
        """Returns the 'keep' and 'tolerates' variants together with run_id, marker_id, biosample_id, variant_id


        :param: variant_tolerate: Boolean: Default False. include "variant_tolerate" variants or not?
        :return: pandas df with columns: run_id, marker_id, biosample_id, variant_id
        """
        # Get portion of variant_known_tsv with either keep or keep+variant_tolerate
        if variant_tolerate:  # get also variant_tolerate variant
            run_marker_biosample_variant_keep_df = self.variant_known_ids_df.loc[
                (self.variant_known_df.action.isin(['keep', 'variant_tolerate']))]
        else:  # do only get keep variants
            run_marker_biosample_variant_keep_df = self.variant_known_ids_df.loc[self.variant_known_df.action == 'keep']
        # run_marker_biosample_variant_keep_df = run_marker_biosample_variant_keep_df.merge(self.variant_known_ids_df,
        #                                         on=['run_id', 'marker_id', 'biosample_id', 'variant_id'])
        # Select run_id, marker_id, biosample_id and variant_id
        run_marker_biosample_variant_keep_df = run_marker_biosample_variant_keep_df[
            ['run_id', 'marker_id', 'biosample_id', 'variant_id']].drop_duplicates(inplace=False)
        # Change variant_id type to int
        run_marker_biosample_variant_keep_df.variant_id = run_marker_biosample_variant_keep_df.variant_id.astype('int')
        return run_marker_biosample_variant_keep_df


    def get_delete_run_marker_biosample_variant_df(self, variant_read_count_df):
        """Returns the 'delete' variants together with run_id, marker_id, biosample_id, variant_id

        :return: pandas df with columns: run_id, marker_id, biosample_id, variant_id
        """

        ##########################################################
        #
        # Get delete variants, that are not keep and/or tolerate in mock or real
        #
        ##########################################################
        # Get mock biosamples
        run_marker_biosample_mock_df = self.variant_known_ids_df.loc[
            self.variant_known_df.biosample_type == 'mock', ['run_id', 'marker_id', 'biosample_id']]
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
        variant_delete_negative_df = self.variant_known_ids_df.loc[self.variant_known_df.biosample_type == 'negative',
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
        variant_delete_real_df = self.variant_known_ids_df.loc[self.variant_known_df.action == 'delete']
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

