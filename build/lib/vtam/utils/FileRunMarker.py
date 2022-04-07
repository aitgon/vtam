import argparse
import pandas
import sqlalchemy

from vtam.models.Marker import Marker
from vtam.models.Run import Run
from vtam.models.SampleInformation import SampleInformation
from vtam.utils.ArgParser import ArgParserChecker
from vtam.utils.NameIdConverter import NameIdConverter


class FileRunMarker(object):

    def __init__(self, tsv_path):
        self.tsv_path = tsv_path

    @staticmethod
    def help():

        return """Input TSV file with two columns and headers 'Run' and 'Marker'.
                                        Example:
                                        Run	Marker
                                        prerun	MFZR
                                        prerun	ZFZR"""

    def get_sample_ids(self, engine):

        sample_lst = []

        record_list = [*self.to_identifier_df(engine).T.to_dict().values()]
        declarative_meta_table = SampleInformation.__table__
        stmt = sqlalchemy.select([declarative_meta_table.c.sample_id]).where(
            declarative_meta_table.c.run_id == sqlalchemy.bindparam('run_id')).where(
            declarative_meta_table.c.marker_id == sqlalchemy.bindparam('marker_id')).distinct()

        with engine.connect() as conn:
            for record in record_list:
                sample_lst = sample_lst + [i[0] for i in conn.execute(stmt, record).fetchall()]

        return [*set(sample_lst)]

    def read_tsv_into_df(self):
        """
        Updated: Mai 10, 2020

        Parameters
        ----------

        Returns
        -------
        pandas.DataFrame

        """

        run_marker_df = pandas.read_csv(self.tsv_path, sep='\t', header=0)
        run_marker_df.columns = run_marker_df.columns.str.lower()
        run_marker_df.rename({'run': 'run_name', 'marker': 'marker_name', 'variant': 'variant_id'}, axis=1, inplace=True)

        return run_marker_df

    def to_identifier_df(self, engine):
        """Takes this file and creates a DF with identifiers.
        Updated: Mai 21, 2020

        Parameters
        ----------
        sqlalchemy engine

        Returns
        -------
        pandas.DataFrame
        DF with ids instead of names and columns run_id, marker_id, sample_id

        """

        run_marker_ids_df = self.read_tsv_into_df().copy()
        run_marker_ids_df['run'] = NameIdConverter(self.read_tsv_into_df().run_name, engine).to_ids(Run)
        run_marker_ids_df['marker'] = NameIdConverter(self.read_tsv_into_df().marker_name, engine).to_ids(Marker)
        run_marker_ids_df.rename({'run': 'run_id', 'marker': 'marker_id'}, axis=1, inplace=True)
        return run_marker_ids_df

    def get_variant_read_count_df(
            self,
            engine,
            variant_read_count_like_model,
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
            stmt_select = sqlalchemy.select(
                [
                    variant_read_count_like_table.c.run_id,
                    variant_read_count_like_table.c.marker_id,
                    variant_read_count_like_table.c.sample_id,
                    variant_read_count_like_table.c.replicate,
                    variant_read_count_like_table.c.variant_id,
                    variant_read_count_like_table.c.read_count]).distinct() .where(
                variant_read_count_like_table.c.run_id == run_id).where(
                variant_read_count_like_table.c.marker_id == marker_id)
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

        return variant_read_count_df

    def check_argument(self):
        """Checks if tsv_path exists, is not empty

        :param tsv_path: Valid non-empty TSV tsv_path tsv_path
        :return: void

        """

        path = ArgParserChecker.check_file_exists_and_is_nonempty(
            self.tsv_path)
        df = pandas.read_csv(path, sep='\t', header=0)
        header = {'Run', 'Marker'}
        if set(df.columns.tolist()) > header:
            raise argparse.ArgumentTypeError(
                "The header of the file {} does not contain these fields: {}!".format(
                    path, header))
        else:
            return path
