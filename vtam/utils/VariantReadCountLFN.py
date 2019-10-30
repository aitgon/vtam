import sys

from vtam.utils.Logger import Logger
from vtam.utils.VTAMexception import VTAMexception


class VariantReadCountLFN(object):
    """
    Takes as input an variant_read_count_df (run_id, marker_id, biosample_id, replicate_id, N_ijk) and returns
    the different LFN calculation, that N_i, N_ik, ...

    N_ijk stands for the read count for each variant_id i, biosample_id j and replicate_id k
    """

    def __init__(self, variant_read_count_df):
        """

        :param variant_read_count_df: DataFrame with columns run_id, marker_id, biosample_id, replicate_id, N_ijk
        """

        try:
            assert sorted(variant_read_count_df.columns.tolist()) == sorted(['run_id', 'marker_id', 'biosample_id', 'replicate_id', 'variant_id', 'read_count'])
        except:
            Logger.instance().error(VTAMexception("This DataFrame is not composed of columns: 'run_id', 'marker_id', "
                                                  "'biosample_id', 'replicate_id', 'variant_id', 'read_count'. The workflow will exit"))
            sys.exit(1)
        self.variant_read_count_df = variant_read_count_df


    def get_N_i_df(self):
        """Returns N_i_df, that is a DataFrame with columns run_id, marker_id, biosample_id, N_ijk
        N_i = sum aggregation of N_ijk over variants i

        """

        N_i_df = self.variant_read_count_df.groupby(by=['run_id', 'marker_id', 'variant_id']).agg({'read_count': sum})\
            .reset_index()
        N_i_df = N_i_df.rename(columns={'read_count': 'N_i'})
        N_i_df.drop_duplicates(inplace=True)

        return N_i_df


    def get_N_ik_df(self):
        """Returns N_ik_df, that is a DataFrame with columns run_id, marker_id, biosample_id, N_ijk
        N_k = sum aggregation of N_ijk over variants i and replicate k

        """

        N_ik_df = self.variant_read_count_df.groupby(by=['run_id', 'marker_id', 'variant_id', 'replicate_id']).agg({'read_count': sum})\
            .reset_index()
        N_ik_df = N_ik_df.rename(columns={'read_count': 'N_ik'})
        N_ik_df.drop_duplicates(inplace=True)

        return N_ik_df


    def get_N_kj_df(self):
        """Returns N_i_df, that is a DataFrame with columns run_id, marker_id, biosample_id, N_ijk
        N_kj = sum aggregation of N_ijk over biosample k and replicate j

        """

        N_kj_df = self.variant_read_count_df.groupby(by=['run_id', 'marker_id', 'biosample_id', 'replicate_id']).agg({'read_count': sum})\
            .reset_index()
        N_kj_df = N_kj_df.rename(columns={'read_count': 'N_kj'})
        N_kj_df.drop_duplicates(inplace=True)

        return N_kj_df

    def filter_singletons(self):
        """Returns a new variant_read_count_df without singletons that is variants i where N_i=1

        :return variant_read_count_df: DataFrame without singletons
        """
        N_i_df = self.get_N_i_df()
        N_i_df = N_i_df.loc[N_i_df.N_i > 1] # get non-singletons
        variant_read_count_df = self.variant_read_count_df.merge(N_i_df, on=['run_id', 'marker_id', 'variant_id'])
        variant_read_count_df.drop('N_i', axis=1, inplace=True)
        # variant_read_count_df.drop_duplicates(inplace=True)
        return variant_read_count_df
