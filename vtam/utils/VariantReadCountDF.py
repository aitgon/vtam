import sys

import pandas
import sqlalchemy

from vtam.models.Biosample import Biosample
from vtam.models.Marker import Marker
from vtam.models.Run import Run
from vtam.utils.Logger import Logger
from vtam.utils.VTAMexception import VTAMexception


class VariantReadCountDF(object):
    """
    Takes as input a variant_read_count_input_df (run_id, marker_id, biosample_id, replicate, N_ijk) and returns
    the different LFN calculation, that is N_i, N_ik, ...

    N_ijk stands for the read count for each variant_id i, biosample_id j and replicate k
    """

    def __init__(self, variant_read_count_df):
        """

        :param variant_read_count_df: DataFrame with columns run_id, marker_id, biosample_id, replicate, N_ijk
        """

        # This is the only order allowed of variant_read_count_input_df columns
        self.column_list = ['run_id', 'marker_id', 'biosample_id', 'replicate', 'variant_id', 'read_count']
        try:
            assert set(variant_read_count_df.columns.tolist()) >= set(self.column_list)
        except:
            Logger.instance().error(VTAMexception("This DataFrame is not composed of columns: 'run_id', 'marker_id', "
                                                  "'biosample_id', 'replicate', 'variant_id', 'read_count'. The workflow will exit"))
            sys.exit(1)

        self.variant_read_count_df = variant_read_count_df

    def get_run_df(self, engine):
        """Based on the SortedReadFile information TSV and Run, returns the run_df

        :return: DataFrame with columns: index, name
        """
        record_list = []
        run_id_list = self.run_id.unique().tolist()
        for run_id in run_id_list:
            stmt_select = sqlalchemy.select([Run.__table__.c.name]).where(Run.__table__.c.id == run_id)
            with engine.connect() as conn:
                for run_name in conn.execute(stmt_select).first():
                    record_list.append({'id': run_id, 'name': run_name})

        run_df = pandas.DataFrame.from_records(record_list, index='id')

        return run_df

    def get_marker_df(self, engine):
        """Returns marker_df for this sample_information_id

        :return: DataFrame with columns: index, name
        """

        record_list = []
        marker_id_list = self.marker_id.unique().tolist()
        for marker_id in marker_id_list:
            stmt_select = sqlalchemy.select([Marker.__table__.c.name]).where(Marker.__table__.c.id == marker_id)
            with engine.connect() as conn:
                for marker_name in conn.execute(stmt_select).first():
                    record_list.append({'id': marker_id, 'name': marker_name})

        marker_df = pandas.DataFrame.from_records(record_list, index='id')

        return marker_df

    def get_biosample_df(self, engine):
        """Returns biosample_df for this sample_information_id

        :return: DataFrame with columns: index, name
        """

        record_list = []
        biosample_id_list = self.variant_read_count_df.biosample_id.unique().tolist()

        for biosample_id in biosample_id_list:
            stmt_select = sqlalchemy.select([Biosample.__table__.c.name]).where(Biosample.__table__.c.id == biosample_id)
            with engine.connect() as conn:
                for biosample_name in conn.execute(stmt_select).first():
                    record_list.append({'id': biosample_id, 'name': biosample_name})

        biosample_df = pandas.DataFrame.from_records(record_list, index='id')

        # Sort biosample names according to fastainfo biosample order
        biosample_id_list = self.to_identifier_df(engine=engine).biosample_id.drop_duplicates(keep='first').tolist()
        biosample_df.index = pandas.Categorical(biosample_df.index, biosample_id_list)
        biosample_df.sort_index(inplace=True)
        biosample_df.index = biosample_df.index.astype(int)

        return biosample_df

    def filter_out_below_global_read_count_threshold(self, global_read_count_threshold):
        """ Returns variants with read_count across all samples with read_count above global_read_count_threshold

        :param global_read_count_threshold: Threshold to get variants. Default throws singletons global_read_count_threshold=2
        :return variant_read_count_input_df: DataFrame without singletons
        """
        N_i_df = self.get_N_i_df()
        N_i_df = N_i_df.loc[N_i_df.N_i >= global_read_count_threshold] # get equal or above global_read_count_threshold
        variant_read_count_df = self.variant_read_count_df.merge(N_i_df, on=['run_id', 'marker_id', 'variant_id'])
        variant_read_count_df.drop('N_i', axis=1, inplace=True)
        # variant_read_count_df.drop_duplicates(inplace=True)
        return variant_read_count_df

    def get_N_i_df(self):
        """Returns N_i_df, that is a DataFrame with columns run_id, marker_id, biosample_id, N_ijk
        N_i = sum aggregation of N_ijk over variants i

        """

        N_i_df = self.variant_read_count_df.groupby(by=['run_id', 'marker_id', 'variant_id']).agg({'read_count': sum})\
            .reset_index()
        N_i_df = N_i_df.rename(columns={'read_count': 'N_i'})
        N_i_df.drop_duplicates(inplace=True)

        return N_i_df

    def get_N_ij_df(self):
        """Returns N_ij_df, that is a DataFrame with columns run_id, marker_id, variant_id, biosample_id, N_ij
        N_ij = sum aggregation of N_ijk over variants i and biosamples j

        """

        N_ij_df = self.variant_read_count_df.groupby(by=['run_id', 'marker_id', 'variant_id', 'biosample_id']).agg({'read_count': sum})\
            .reset_index()
        N_ij_df = N_ij_df.rename(columns={'read_count': 'N_ij'})
        N_ij_df.drop_duplicates(inplace=True)

        return N_ij_df

    def get_N_ik_df(self):
        """Returns N_ik_df, that is a DataFrame with columns run_id, marker_id, biosample_id, N_ijk
        N_k = sum aggregation of N_ijk over variants i and replicate k

        """

        N_ik_df = self.variant_read_count_df.groupby(by=['run_id', 'marker_id', 'variant_id', 'replicate']).agg({'read_count': sum})\
            .reset_index()
        N_ik_df = N_ik_df.rename(columns={'read_count': 'N_ik'})
        N_ik_df.drop_duplicates(inplace=True)

        return N_ik_df

    def get_N_jk_df(self):

        """Returns N_i_df, that is a DataFrame with columns run_id, marker_id, biosample_id, N_ijk
        N_kj = sum aggregation of N_ijk over biosample j and replicate k

        """

        N_jk_df = self.variant_read_count_df.groupby(by=['run_id', 'marker_id', 'biosample_id', 'replicate']).agg({'read_count': sum})\
            .reset_index()
        N_jk_df = N_jk_df.rename(columns={'read_count': 'N_jk'})
        N_jk_df.drop_duplicates(inplace=True)

        return N_jk_df

    def convert_to_record_list(self):

        """Convert DF to list of dictionaries to use in an sqlalchemy core insert"""

        records = []
        # import pdb; pdb.set_trace()
        for row in self.variant_read_count_df.itertuples():
            run_id = row.run_id
            marker_id = row.marker_id
            biosample_id = row.biosample_id
            variant_id = row.variant_id
            read_count = row.read_count
            instance = {'run_id': run_id, 'marker_id': marker_id,
                                                         'variant_id': variant_id, 'biosample_id': biosample_id,
                                                         'read_count': read_count}
            if 'filter_delete' in dir(row):
                instance['filter_delete'] = row.filter_delete
            if 'filter_id' in dir(row):
                instance['filter_id'] = row.filter_id
            if 'replicate' in dir(row):
                instance['replicate'] = row.replicate
            if 'replicate_count' in dir(row):
                instance['replicate_count'] = row.replicate_count
            if 'read_count_average' in dir(row):
                instance['read_count_average'] = row.read_count_average
            records.append(instance)
        return records
