import argparse

import pandas
import sqlalchemy

from vtam.models.Marker import Marker
from vtam.models.Run import Run
from vtam.utils.ArgParser import ArgParserChecker


class RunMarkerFile(object):

    def __init__(self, tsv_path):
        self.tsv_path = tsv_path

    @staticmethod
    def help():

        return """Input TSV file with two columns and headers 'Run' and 'Marker'.
                                        Example:
                                        Run	Marker
                                        prerun	MFZR
                                        prerun	ZFZR"""

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

        return run_marker_df

    def to_identifier_df(self, engine):
        """Takes this file and creates a DF with identifiers

        Parameters
        ----------
        engine sqlalchemy engine

        Returns
        -------
        pandas.DataFrame
        DF with ids instead of names and columns run_id, marker_id, biosample_id

        """

        sample_info_df = pandas.DataFrame()
        sample_info_df['run_id'] = None
        sample_info_df['marker_id'] = None
        with engine.connect() as conn:
            for k, row in self.read_tsv_into_df().iterrows():
                marker_name = row.marker
                run_name = row.run
                # get run_id ###########
                stmt_select_run_id = sqlalchemy.select([Run.__table__.c.id]).where(Run.__table__.c.name == run_name)
                run_id = conn.execute(stmt_select_run_id).first()[0]
                # get marker_id ###########
                stmt_select_marker_id = sqlalchemy.select([Marker.__table__.c.id]).where(Marker.__table__.c.name == marker_name)
                marker_id = conn.execute(stmt_select_marker_id).first()[0]
                new_row_dic = dict(row)
                new_row_dic['run_id'] = run_id
                new_row_dic['marker_id'] = marker_id
                del new_row_dic["run"]
                del new_row_dic["marker"]

                sample_info_df = sample_info_df.append(pandas.DataFrame(new_row_dic, index=[0]), sort=False)

        sample_info_df.columns = sample_info_df.columns.str.lower()
        return sample_info_df

    def get_variant_read_count_df(self, engine, variant_read_count_like_model, filter_id=None):
        """Based on the SortedReadFile samples and the variant_read_count_model, returns the variant_read_count_input_df

        :param variant_read_count_like_model: SQLalchemy models with columns: run_id, marker_id, biosample_id, replicate, variant_id, read_count
        :param filter_id:
        :return: DataFrame with columns: run_id, marker_id, biosample_id, replicate, variant_id, read_count
        """

        variant_read_count_like_table = variant_read_count_like_model.__table__

        variant_read_count_list = []
        # for sample_instance in self.get_fasta_information_record_list():
        for sample_instance_row in self.to_identifier_df(engine=engine).itertuples():
            run_id = sample_instance_row.run_id
            marker_id = sample_instance_row.marker_id
            stmt_select = sqlalchemy.select([
                variant_read_count_like_table.c.run_id, variant_read_count_like_table.c.marker_id,
                variant_read_count_like_table.c.biosample_id, variant_read_count_like_table.c.replicate,
                variant_read_count_like_table.c.variant_id, variant_read_count_like_table.c.read_count]).distinct()\
                .where(variant_read_count_like_table.c.run_id == run_id).where(variant_read_count_like_table.c.marker_id == marker_id)
            # Used for filters tables where filter_delete attribute exists
            if 'filter_delete' in [column.key for column in variant_read_count_like_table.columns]:
                stmt_select = stmt_select.where(variant_read_count_like_table.c.filter_delete == 0)
            # used for filter lfn where filter_id = 8 is necessary (do not pass all filters)
            if not filter_id is None:
                stmt_select = stmt_select.where(variant_read_count_like_table.c.filter_id == filter_id)
            with engine.connect() as conn2:
                for row in conn2.execute(stmt_select).fetchall():
                    variant_read_count_list.append(row)
        #
        variant_read_count_df = pandas.DataFrame.from_records(variant_read_count_list,
            columns=['run_id', 'marker_id', 'biosample_id', 'replicate', 'variant_id', 'read_count'])

        return variant_read_count_df

    def check_argument(self):

        """Checks if tsv_path exists, is not empty

        :param tsv_path: Valid non-empty TSV tsv_path tsv_path
        :return: void

        """

        path = ArgParserChecker.check_file_exists_and_is_nonempty(self.tsv_path)
        df = pandas.read_csv(path, sep='\t', header=0)
        header = {'Run', 'Marker'}
        if set(df.columns.tolist()) > header:
            raise argparse.ArgumentTypeError("The header of the file {} does not contain these fields: {}!".format(path, header))
        else:
            return path

