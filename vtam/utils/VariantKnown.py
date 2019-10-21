import os
import pandas
import sys

import sqlalchemy
from sqlalchemy import select

from vtam import Logger, VTAMexception
from vtam.utils.FastaInfo import FastaInfo


class VariantKnown(object):

    def __init__(self, variant_known_tsv, fasta_info_tsv, engine, variant_model, run_model, marker_model, biosample_model, replicate_model):
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
        self.replicate_model = replicate_model
        self.fasta_info_tsv = fasta_info_tsv

        ################################################################################################################
        #
        # Create a DF
        #
        ################################################################################################################
        self.variant_known_df = pandas.read_csv(self.variant_known_tsv, sep="\t", header=0, \
                                              names=['marker_name', 'run_name', 'biosample_name', 'biosample_type',
                                                     'variant_id', 'action', 'variant_sequence', 'note'], index_col=False)
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
        """Returns a list of dictionnaries with run_id, marker_id, biosample_id and replicate_id entries (See return)

        :return: list of dictionnaries: [{'run_id': 1, 'marker_id': 1, 'biosample_id': 1, 'replicate_id': 1}, {'run_id': 1, ...
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
                # get replicate_id ###########
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
                    Logger.instance().error(VTAMexception("Error: Verify the IDs and sequence of the variants in the --variant_known file"))
                    sys.exit(1)



    def __are_known_variants_coherent_with_fasta_info_file(self):

        fasta_info = FastaInfo(fasta_info_tsv=self.fasta_info_tsv, engine=self.engine)
        fasta_info_records = fasta_info.get_ids_of_run_marker_biosample_replicate(self.engine, self.run_model,
                                                          self.marker_model, self.biosample_model, self.replicate_model)
        fasta_info_df = pandas.DataFrame.from_records(data=fasta_info_records)

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
                assert (fasta_info_df.loc[
                    (fasta_info_df['biosample_id'] == biosample_id) & (fasta_info_df['marker_id'] == marker_id) & (
                                fasta_info_df['run_id'] == run_id)]).shape[0] > 0
            except AssertionError:
                Logger.instance().error(VTAMexception("Error: Verify in the --variant_known file that run_id, marker_id and biosample_id"
                                                      "are defined in the --fasta_info file"))
                sys.exit(1)
