import sys

from vtam.utils.Logger import Logger
from vtam.utils.VTAMexception import VTAMexception


class DataframeVariantReadCountLike(object):
    """
    Takes as input a variant_read_count_input_df (run_id, marker_id, sample_id, replicate, N_ijk) and returns
    the different LFN calculation, that is N_i, N_ik, ...

    N_ijk stands for the read count for each variant_id i, sample_id j and replicate k
    """

    def __init__(self, variant_read_count_df):
        """

        :param variant_read_count_df: DataFrame with columns run_id, marker_id, sample_id, replicate, N_ijk
        """

        # This is the only order allowed of variant_read_count_input_df columns
        self.column_list = [
            'run_id',
            'marker_id',
            'sample_id',
            'replicate',
            'variant_id',
            'read_count']
        try:
            assert set(
                variant_read_count_df.columns.tolist()) >= set(
                self.column_list)
        except BaseException:
            Logger.instance().error(
                VTAMexception(
                    "This DataFrame is not composed of columns: 'run_id', 'marker_id', "
                    "'sample_id', 'replicate', 'variant_id', 'read_count'. The workflow will exit"))
            sys.exit(1)

        self.variant_read_count_df = variant_read_count_df

    def filter_out_below_global_read_count_cutoff(
            self, global_read_count_cutoff):
        """ Returns variants with read_count across all samples with read_count above global_read_count_cutoff

        :param global_read_count_cutoff: Threshold to get variants. Default throws singletons global_read_count_cutoff=2
        :return variant_read_count_input_df: DataFrame without singletons
        """
        N_i_df = self.get_N_i_df()
        #  get equal or above global_read_count_cutoff
        N_i_df = N_i_df.loc[N_i_df.N_i >= global_read_count_cutoff]
        variant_read_count_df = self.variant_read_count_df.merge(
            N_i_df, on=['run_id', 'marker_id', 'variant_id'])
        variant_read_count_df.drop('N_i', axis=1, inplace=True)
        # nijk_df.drop_duplicates(inplace=True)
        return variant_read_count_df

    def get_N_i_df(self):
        """Returns N_i_df, that is a DataFrame with columns run_id, marker_id, sample_id, N_ijk
        N_i = sum aggregation of N_ijk over variants i

        """

        N_i_df = self.variant_read_count_df.groupby(
            by=['run_id', 'marker_id', 'variant_id']).agg({'read_count': sum}) .reset_index()
        N_i_df = N_i_df.rename(columns={'read_count': 'N_i'})
        N_i_df.drop_duplicates(inplace=True)

        return N_i_df

    def get_N_ij_df(self):
        """Returns N_ij_df, that is a DataFrame with columns run_id, marker_id, variant_id, sample_id, N_ij
        N_ij = sum aggregation of N_ijk over variants i and samples j

        """

        N_ij_df = self.variant_read_count_df.groupby(
            by=['run_id', 'marker_id', 'variant_id', 'sample_id']).agg({'read_count': sum}) .reset_index()
        N_ij_df = N_ij_df.rename(columns={'read_count': 'N_ij'})
        N_ij_df.drop_duplicates(inplace=True)

        return N_ij_df

    def get_N_ik_df(self):
        """Returns N_ik_df, that is a DataFrame with columns run_id, marker_id, sample_id, N_ijk
        N_k = sum aggregation of N_ijk over variants i and replicate k

        """

        N_ik_df = self.variant_read_count_df.groupby(
            by=['run_id', 'marker_id', 'variant_id', 'replicate']).agg({'read_count': sum}) .reset_index()
        N_ik_df = N_ik_df.rename(columns={'read_count': 'N_ik'})
        N_ik_df.drop_duplicates(inplace=True)

        return N_ik_df

    def get_N_jk_df(self):
        """Returns N_i_df, that is a DataFrame with columns run_id, marker_id, sample_id, N_ijk
        N_kj = sum aggregation of N_ijk over sample j and replicate k

        """

        N_jk_df = self.variant_read_count_df.groupby(
            by=['run_id', 'marker_id', 'sample_id', 'replicate']).agg({'read_count': sum}) .reset_index()
        N_jk_df = N_jk_df.rename(columns={'read_count': 'N_jk'})
        N_jk_df.drop_duplicates(inplace=True)

        return N_jk_df

    def to_sql(self, engine, variant_read_count_like_model):
        """Convert DF to list of dictionaries to use in an sqlalchemy core insert"""

        record_list = []

        for row in self.variant_read_count_df.itertuples():

            run_id = row.run_id
            marker_id = row.marker_id
            sample_id = row.sample_id
            variant_id = row.variant_id
            read_count = row.read_count
            instance = {'run_id': run_id, 'marker_id': marker_id,
                        'variant_id': variant_id, 'sample_id': sample_id,
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
            record_list.append(instance)

        with engine.connect() as conn:

            # Insert new instances
            conn.execute(
                variant_read_count_like_model.__table__.insert(),
                record_list)
