from vtam.utils.RunnerFilterLFN import RunnerFilterLFN
from vtam.utils.RunnerFilterMinReplicateNumber import RunnerFilterMinReplicateNumber


class RunnerFilterLFNreplicateRemain:
    """This is the filter(s) used by OptimizeLFNreadCountAndLFNvariant"""

    def __init__(self, nijk_df, lfn_ni_cutoff, lfn_nik_cutoff, lfn_njk_cutoff, lfn_nijk_cutoff,
                           min_replicate_number):

        self.nijk_df = nijk_df
        self.lfn_ni_cutoff = lfn_ni_cutoff
        self.lfn_nik_cutoff = lfn_nik_cutoff
        self.lfn_njk_cutoff = lfn_njk_cutoff
        self.lfn_nijk_cutoff = lfn_nijk_cutoff
        self.min_replicate_number = min_replicate_number

    def count_keep_delete(self, known_occurrences_df):

        keep_run_marker_sample_variant_df = known_occurrences_df.loc[
            known_occurrences_df.action == 'keep', ['run_id', 'marker_id', 'sample_id',
                                                           'variant_id']]

        delete_run_marker_sample_variant_df = known_occurrences_df.loc[
            known_occurrences_df.action == 'delete', ['run_id', 'marker_id', 'sample_id',
                                                           'variant_id']]

        nijk_remain_df = self.get_nijk_remain_df()
        count_keep = nijk_remain_df.merge(keep_run_marker_sample_variant_df,
                                          on=['run_id', 'marker_id', 'sample_id',
                                              'variant_id']).drop_duplicates(inplace=False).shape[0]
        count_delete = nijk_remain_df.merge(
            delete_run_marker_sample_variant_df,
            on=['run_id', 'marker_id', 'sample_id', 'variant_id']).drop_duplicates(
            inplace=False).shape[0]

        return count_keep, count_delete

    def get_nijk_remain_df(self):

        lfn_filter_runner = RunnerFilterLFN(self.nijk_df)

        ############################################################################################
        #
        # Filter lfn_variant nijk_ni
        #
        ############################################################################################

        if self.lfn_nik_cutoff is None:  # optimize lfn variant

            lfn_filter_runner.mark_delete_lfn_per_Ni_or_Nik_or_Njk(lfn_denominator='N_i', cutoff=self.lfn_ni_cutoff)

        else:  # optimize lfn variant replicate

            lfn_filter_runner.mark_delete_lfn_per_Ni_or_Nik_or_Njk(lfn_denominator='N_ik', cutoff=self.lfn_nik_cutoff)

        ############################################################################################
        #
        # Filter lfn_sample_replicate
        #
        ############################################################################################

        lfn_filter_runner.mark_delete_lfn_per_Ni_or_Nik_or_Njk(lfn_denominator='N_jk',
                                                               cutoff=self.lfn_njk_cutoff)

        ############################################################################################
        #
        # Filter absolute read count
        #
        ############################################################################################

        lfn_filter_runner.mark_delete_lfn_absolute_read_count(
            lfn_read_count_cutoff=self.lfn_nijk_cutoff)

        ############################################################################################
        #
        # mark_delete_lfn_do_not_pass_all_filters
        #
        ############################################################################################

        lfn_filter_runner.mark_delete_lfn_do_not_pass_all_filters()

        nijk_remain_df = lfn_filter_runner.variant_read_count_filter_delete_df

        nijk_remain_df = nijk_remain_df.loc[
            (nijk_remain_df.filter_id == 8) &
            (nijk_remain_df.filter_delete == 0)]

        del (lfn_filter_runner)

        ############################################################################################
        #
        # RunnerFilterMinReplicateNumber
        #
        ############################################################################################

        nijk_remain_df = RunnerFilterMinReplicateNumber(
            nijk_remain_df).get_variant_read_count_delete_df(self.min_replicate_number)
        nijk_remain_df = nijk_remain_df.loc[
            (nijk_remain_df.filter_delete == 0)]
        nijk_remain_df.drop('filter_delete', axis=1, inplace=True)

        # Delete object

        return nijk_remain_df[
            ['run_id', 'marker_id', 'sample_id', 'variant_id']].copy().drop_duplicates(
            inplace=False)

