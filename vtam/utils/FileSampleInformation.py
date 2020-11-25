import argparse
import os
import sys

import pandas
import sqlalchemy

from vtam.utils.Logger import Logger
from vtam.utils.VTAMexception import VTAMexception

from vtam.models.Run import Run
from vtam.models.Marker import Marker
from vtam.models.Sample import Sample
from vtam.models.Variant import Variant
from vtam.models.SampleInformation import SampleInformation
from vtam.models.SortedReadFile import SortedReadFile


class FileSampleInformation:
    """Sample information file (paired fastq, merged fasta and sorted read fasta)"""

    def __init__(self, tsv_path):
        self.tsv_path = tsv_path

    def delete_from_db(self, engine, variant_read_count_like_model):
        """Convert DF to list of dictionaries to use in an sqlalchemy core insert"""

        sample_information_df = self.to_identifier_df(engine=engine)
        sample_record_list = sample_information_df.to_dict('records')
        sample_column_list = sample_information_df.columns.tolist()

        with engine.connect() as conn:
            stmt = variant_read_count_like_model.__table__.delete()
            stmt = stmt.where(
                variant_read_count_like_model.__table__.c.run_id == sqlalchemy.bindparam('run_id'))
            stmt = stmt.where(
                variant_read_count_like_model.__table__.c.marker_id == sqlalchemy.bindparam('marker_id'))
            if 'sample_id' in sample_column_list:
                stmt = stmt.where(
                    variant_read_count_like_model.__table__.c.sample_id == sqlalchemy.bindparam('sample_id'))
            if 'replicate' in sample_column_list and 'replicate' in [
                    col.key for col in variant_read_count_like_model.__table__.columns]:
                stmt = stmt.where(
                    variant_read_count_like_model.__table__.c.replicate == sqlalchemy.bindparam('replicate'))
            conn.execute(stmt, sample_record_list)

    def get_nijk_df(
            self,
            variant_read_count_like_model,
            engine,
            filter_id=None):
        """Based on the SortedReadFile samples and the variant_read_count_model, returns the variant_read_count_input_df

        :param variant_read_count_like_model: SQLalchemy models with columns: run_id, marker_id, sample_id, replicate, variant_id, read_count
        :param filter_id:
        :return: DataFrame with columns: run_id, marker_id, sample_id, replicate, variant_id, read_count
        """

        variant_read_count_like_table = variant_read_count_like_model.__table__

        variant_read_count_list = []
        # for sample_instance in self.get_fasta_information_record_list():
        for sample_instance_row in self.to_identifier_df(
                engine=engine).itertuples():
            run_id = sample_instance_row.run_id
            marker_id = sample_instance_row.marker_id
            sample_id = sample_instance_row.sample_id
            replicate = sample_instance_row.replicate
            stmt_select = sqlalchemy.select(
                [
                    variant_read_count_like_table.c.run_id,
                    variant_read_count_like_table.c.marker_id,
                    variant_read_count_like_table.c.sample_id,
                    variant_read_count_like_table.c.replicate,
                    variant_read_count_like_table.c.variant_id,
                    variant_read_count_like_table.c.read_count]).distinct() .where(
                variant_read_count_like_table.c.run_id == run_id).where(
                variant_read_count_like_table.c.marker_id == marker_id) .where(
                    variant_read_count_like_table.c.sample_id == sample_id).where(
                        variant_read_count_like_table.c.replicate == replicate)
            # Used for filters tables where filter_delete attribute exists
            if 'filter_delete' in [
                    column.key for column in variant_read_count_like_table.columns]:
                stmt_select = stmt_select.where(
                    variant_read_count_like_table.c.filter_delete == 0)
            # used for filter lfn where filter_id = 8 is necessary (do not pass
            # all filters)
            if filter_id is not None:
                stmt_select = stmt_select.where(
                    variant_read_count_like_table.c.filter_id == filter_id)
            with engine.connect() as conn2:
                for row in conn2.execute(stmt_select).fetchall():
                    variant_read_count_list.append(row)
        #
        variant_read_count_df = pandas.DataFrame.from_records(
            variant_read_count_list,
            columns=[
                'run_id',
                'marker_id',
                'sample_id',
                'replicate',
                'variant_id',
                'read_count'])

        # Exit if no variants for analysis
        try:
            assert variant_read_count_df.shape[0] > 0
        except AssertionError:
            Logger.instance().warning(
                VTAMexception(
                    "No variants available after this filter. "
                    "The pipeline will stop here.".format(
                        self.__class__.__name__)))
            sys.exit(0)
        return variant_read_count_df

    def to_identifier_df(self, engine):
        """Takes the sample information stats_df and replaces names with ids

        Returns
        -------
        pandas.DataFrame
        DF with ids instead of names and columns run_id, marker_id, sample_id

        """

        sample_info_df = pandas.DataFrame()
        sample_info_df['run_id'] = None
        sample_info_df['marker_id'] = None
        sample_info_df['sample_id'] = None
        with engine.connect() as conn:
            for k, row in self.read_tsv_into_df().iterrows():
                marker_name = row.marker
                run_name = row.run
                sample_name = row['sample']  # fixed confusion of sample method and variable
                replicate = int(row.replicate)
                # get run_id ###########
                stmt_select_run_id = sqlalchemy.select(
                    [Run.__table__.c.id]).where(Run.__table__.c.name == run_name)
                run_id = conn.execute(stmt_select_run_id).first()[0]
                # get marker_id ###########
                stmt_select_marker_id = sqlalchemy.select([Marker.__table__.c.id]).where(
                    Marker.__table__.c.name == marker_name)
                marker_id = conn.execute(stmt_select_marker_id).first()[0]
                # get sample_id ###########
                stmt_select_sample_id = sqlalchemy.select([Sample.__table__.c.id]).where(
                    Sample.__table__.c.name == sample_name)
                sample_id = conn.execute(
                    stmt_select_sample_id).first()[0]
                new_row_dic = dict(row)
                new_row_dic['run_id'] = run_id
                new_row_dic['marker_id'] = marker_id
                new_row_dic['sample_id'] = sample_id
                new_row_dic['replicate'] = replicate
                del new_row_dic["run"]
                del new_row_dic["marker"]
                del new_row_dic["sample"]

                sample_info_df = sample_info_df.append(
                    pandas.DataFrame(new_row_dic, index=[0]), sort=False)

        sample_info_df.columns = sample_info_df.columns.str.lower()
        return sample_info_df

    def read_tsv_into_df(self):
        """
        Updated: Mai 10, 2020

        Parameters
        ----------
        tsv_path: str
            Path to the information files in TSV format

        Returns
        -------
        pandas.DataFrame

        """

        fastqinfo_df = pandas.read_csv(self.tsv_path, sep='\t', header=0)
        fastqinfo_df.columns = fastqinfo_df.columns.str.lower()

        return fastqinfo_df

    def check_args(self, header):
        """
        Checks that the this input file is in the correct format. It is also used by the argparser.
        Updated: Oct 6, 2020

        Parameters
        ----------
        header: list
            List of header fields in lower case of FileSampleInformation

        Returns
        -------

        """

        """Check tsv_path format

        :param tsv_path: Valid non-empty file tsv_path
        :return: void

        """

        if not os.path.isfile(self.tsv_path):
            raise argparse.ArgumentTypeError(
                "The file '{}' does not exist. Please fix it.".format(
                    self.tsv_path))
        elif not os.stat(self.tsv_path).st_size > 0:
            raise argparse.ArgumentTypeError(
                "The file '{}' is empty!".format(self.tsv_path))

        sample_info_df = self.read_tsv_into_df()

        # Verify that the file contains at least the 'header_lower' columns
        if set(sample_info_df.columns) >= set(header):
            # Verify if combination of run, marker, sample, replicates values are unique
            if sample_info_df[['run', 'marker', 'sample', 'replicate']].drop_duplicates().shape[0] == sample_info_df[['run', 'marker', 'sample', 'replicate']].shape[0]:
                return self.tsv_path  # return the tsv_path
            else:
                raise argparse.ArgumentTypeError(
                    "The combinations of the columns 'run, marker, sample, replicate' must be unique in this file '{}'."
                    " Please fix it.".format(self.tsv_path))
        else:
            raise argparse.ArgumentTypeError(
                "The file '{}' must include at least these columns: {}. Please fix it.".format(
                    self.tsv_path, header))

    def to_sqlite(self, session):
        """Takes the sample information stats_df and replaces names with ids

        Returns
        -------
        pandas.DataFrame
        DF with ids instead of names and columns run_id, marker_id, sample_id

        """

        for row in self.read_tsv_into_df().itertuples():
            run_name = row.run
            marker_name = row.marker
            sample_name = row.sample
            replicate = row.replicate
            sorted_read_file = row.sortedfasta
            #
            # Insert run_name
            run_obj = {'name': run_name}
            run_instance = FileSampleInformation.get_or_create(
                session, Run, **run_obj)
            run_id = run_instance.id
            #
            # Insert marker_id
            marker_obj = {'name': marker_name}
            marker_instance = FileSampleInformation.get_or_create(
                session, Marker, **marker_obj)
            marker_id = marker_instance.id
            #
            # Insert Biosamples
            sample_obj = {'name': sample_name}
            sample_instance = FileSampleInformation.get_or_create(
                session, Sample, **sample_obj)
            sample_id = sample_instance.id
            #
            # Insert file output
            fasta_obj = {'name': sorted_read_file, 'run_id': run_id}
            fasta_instance = FileSampleInformation.get_or_create(
                session, SortedReadFile, **fasta_obj)
            fasta_id = fasta_instance.id

            # Insert sample_information
            sample_information_obj = {
                'sample_id': sample_id,
                'replicate': replicate,
                'run_id': run_id}
            sample_information_obj['sortedreadfile_id'] = fasta_id
            sample_information_obj['marker_id'] = marker_id
            FileSampleInformation.get_or_create(
                session, SampleInformation, **sample_information_obj)

    def get_variant_df(
            self,
            variant_read_count_like_model,
            engine,
            filter_id=None):
        """Based on the SortedReadFile information TSV and variant_model, returns the variant_df

        :return: DataFrame with columns: index, sequence
        """

        variant_read_count_like_df = self.get_nijk_df(
            variant_read_count_like_model, engine, filter_id)

        record_list = []
        variant_model_table = Variant.__table__
        for variant_id in variant_read_count_like_df.variant_id.unique().tolist():
            with engine.connect() as conn:
                stmt_select = sqlalchemy.select([variant_model_table.c.sequence]).where(
                    variant_model_table.c.id == variant_id)
                for sequence in conn.execute(stmt_select).first():
                    record_list.append(
                        {'id': variant_id, 'sequence': sequence})

        variant_df = pandas.DataFrame.from_records(record_list, index='id')

        return variant_df

    @staticmethod
    def get_or_create(session, model, **kwargs):
        instance = session.query(model).filter_by(**kwargs).first()
        if instance:  # get
            return instance
        else:  # create
            instance = model(**kwargs)
            session.add(instance)
            session.commit()
            return instance
