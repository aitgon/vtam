import math
import multiprocessing
import os
import pandas
import pathlib

from vtam.utils.PathManager import PathManager
from vtam.utils.RunnerVSearch import RunnerVSearch
from vtam.utils.DataframeVariant import DataframeVariant
from vtam.utils.DataframeVariantReadCountLike import DataframeVariantReadCountLike


class RunnerFilterPCRerror(object):
    """Has attributes and methods to run_name the PCR error Filter"""

    def __init__(
            self,
            variant_expected_df,
            variant_unexpected_df,
            variant_read_count_df):
        """
        Initiates object for the PCR error filter

        :param variant_expected_df: DataFrame (id, sequence) with expected variants
        :param variant_unexpected_df: DataFrame (id, sequence) with unexpected variants
        :param variant_read_count_df: DataFrame (run_id, marker_id, sample_id, replicate, variant_id, read_count)
        """
        self.__variant_expected_df = variant_expected_df
        self.__variant_unexpected_df = variant_unexpected_df
        self.__variant_read_count_df = variant_read_count_df
        self.__tmp_dir = os.path.join(PathManager.instance().get_tempdir(), self.__class__.__name__)
        pathlib.Path(self.__tmp_dir).mkdir(parents=True, exist_ok=True)

    def get_variant_read_count_delete_df(self, pcr_error_var_prop):

        variant_unexpected_to_expected_ratio_df = self.get_variant_unexpected_to_expected_ratio_df()

        # Initiates filter_output_df
        filter_output_df = self.__variant_read_count_df.copy()
        filter_output_df['filter_delete'] = False

        for row in variant_unexpected_to_expected_ratio_df.itertuples():
            if float(getattr(row, 'N_ij_unexpected_to_expected_ratio')
                     ) < pcr_error_var_prop:
                filter_output_df.loc[(filter_output_df['run_id'] == row.run_id)
                                     & (filter_output_df['marker_id'] == row.marker_id)
                                     & (filter_output_df['sample_id'] == row.sample_id)
                                     & (filter_output_df['variant_id'] == row.variant_id_unexpected), 'filter_delete'] = True
        return filter_output_df

    def get_vsearch_alignement_df(self):
        """
        This function runs vsearch to detect PCR errors (1 mism or gap) between the "db" and the "query" sets

        Returns: Pandas DataFrame with output of vsearch and these columnts: query, target, alnlen, ids, mism, gaps
        """

        # length of smallest sequence
        length_min = min(
            self.__variant_expected_df.sequence.apply(len).tolist() +
            self.__variant_unexpected_df.sequence.apply(len).tolist())
        # calcul identity
        identity = math.floor((length_min - 1) / length_min * 100) / 100

        #
        ###################################################################
        # 5-1. Make a fasta_path file with all variants of the sample or replicate
        ###################################################################

        variant_expected_fasta_path = os.path.join(
            self.__tmp_dir, '{}.fasta'.format("variant_expected"))
        variant_expected_df_utils_obj = DataframeVariant(
            variant_df=self.__variant_expected_df)
        variant_expected_df_utils_obj.to_fasta(
            fasta_path=variant_expected_fasta_path)

        variant_unexpected_fasta_path = os.path.join(
            self.__tmp_dir, '{}.fasta'.format("variant_unexpected"))
        variant_unexpected_df_utils_obj = DataframeVariant(
            variant_df=self.__variant_unexpected_df)
        variant_unexpected_df_utils_obj.to_fasta(
            fasta_path=variant_unexpected_fasta_path)

        #
        # Create object and run_name vsearch
        if os.getenv('VTAM_THREADS') is None:
            num_threads = multiprocessing.cpu_count()
        else:
            num_threads = int(os.getenv('VTAM_THREADS'))
        vsearch_pcr_error_tsv = os.path.join(
            self.__tmp_dir, '{}.tsv'.format("vsearch_pcr_error"))
        vsearch_parameters = {
            'db': variant_expected_fasta_path,
            'usearch_global': variant_unexpected_fasta_path,
            'id': str(identity),
            'maxrejects': 0,
            'maxaccepts': 0,
            'userout': vsearch_pcr_error_tsv,
            'userfields': "query+target+alnlen+ids+mism+gaps",
            'threads': num_threads,
        }
        vsearch_cluster = RunnerVSearch(parameters=vsearch_parameters)
        vsearch_cluster.run()

        column_names = [
            'variant_id_unexpected',
            'variant_id_expected',
            'alnlen',
            'ids',
            'mism',
            'gaps']
        vsearch_alignement_df = pandas.read_csv(
            vsearch_pcr_error_tsv, sep='\t', names=column_names)
        return vsearch_alignement_df

    def get_variant_unexpected_to_expected_ratio_df(self):
        """Creates a DF with these columns
        ['run_id', 'marker_id', 'sample_id', 'variant_id_expected', 'N_ij_expected', 'variant_id_unexpected',
        'N_ij_unexpected', 'N_ij_unexpected_to_expected_ratio']

        """

        ############################################################################################
        #
        # Get variant pairs with 1 sequence difference (mism and/or gaps)
        #
        #############################################################################################

        pcr_error_df = self.get_vsearch_alignement_df()
        # Add up mismatch and gap
        pcr_error_df[
            'sum_mism_gaps'] = pcr_error_df.mism + pcr_error_df.gaps
        pcr_error_df = pcr_error_df.loc[
            pcr_error_df.sum_mism_gaps == 1, ['variant_id_expected', 'variant_id_unexpected']]

        ############################################################################################
        #
        # Append N_ij_expected and N_ij_unexpected
        #
        ############################################################################################

        variant_read_count_lfn_instance = DataframeVariantReadCountLike(self.__variant_read_count_df)
        N_ij_df = variant_read_count_lfn_instance.get_N_ij_df()

        pcr_error_df = pcr_error_df.merge(N_ij_df, left_on=['variant_id_expected'], right_on=['variant_id'])
        pcr_error_df.rename(columns={'N_ij': 'N_ij_expected'}, inplace=True)
        pcr_error_df.drop('variant_id', axis=1, inplace=True)

        pcr_error_df = pcr_error_df.merge(N_ij_df, left_on=[
            'run_id', 'marker_id', 'sample_id', 'variant_id_unexpected'], right_on=['run_id', 'marker_id', 'sample_id', 'variant_id'])
        pcr_error_df.rename(columns={'N_ij': 'N_ij_unexpected'}, inplace=True)
        pcr_error_df.drop('variant_id', axis=1, inplace=True)

        ############################################################################################
        #
        # Ratio variant_unexpected_to_expected_ratio_df
        #
        ############################################################################################

        # Add two column for the two expected ratio cases ratio 1 and ratio 2
        pcr_error_df['N_ij_unexpected_to_expected_ratio'] = pcr_error_df['N_ij_unexpected'] \
                                                            / pcr_error_df['N_ij_expected']
        # reorder columns
        pcr_error_df = pcr_error_df[[
            'run_id', 'marker_id', 'sample_id', 'variant_id_expected', 'N_ij_expected',
            'variant_id_unexpected', 'N_ij_unexpected', 'N_ij_unexpected_to_expected_ratio']]

        return pcr_error_df
