from vtam.utils.FilterLFNRunner import FilterLFNrunner
from vtam.utils.FilterMinReplicateNumberRunner import FilterMinReplicateNumberRunner


class FilterLFNreplicateRemainRunner:
    """This is the filter(s) used by OptimizeLFNreadCountAndLFNvariant"""

    def __init__(self, nijk_df):

        self.nijk_df = nijk_df

    def get_nijk_remain_df(self, lfn_ni_cutoff, lfn_nik_cutoff, lfn_njk_cutoff, lfn_nijk_cutoff,
                           min_replicate_number):

        lfn_filter_runner = FilterLFNrunner(self.nijk_df)

        ############################################################################################
        #
        # Filter lfn_variant nijk_ni
        #
        ############################################################################################

        if lfn_nik_cutoff is None:  # optimize lfn variant

            lfn_filter_runner.mark_delete_lfn_per_Ni_or_Nik_or_Njk(lfn_denominator='N_i', threshold=lfn_ni_cutoff)

        else:  # optimize lfn variant replicate

            lfn_filter_runner.mark_delete_lfn_per_Ni_or_Nik_or_Njk(lfn_denominator='N_ik', threshold=lfn_nik_cutoff)

        ############################################################################################
        #
        # Filter lfn_biosample_replicate
        #
        ############################################################################################

        lfn_filter_runner.mark_delete_lfn_per_Ni_or_Nik_or_Njk(lfn_denominator='N_jk',
                                                               threshold=lfn_njk_cutoff)

        ############################################################################################
        #
        # Filter absolute read count
        #
        ############################################################################################

        lfn_filter_runner.mark_delete_lfn_absolute_read_count(
            lfn_read_count_threshold=lfn_nijk_cutoff)

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
        # FilterMinReplicateNumberRunner
        #
        ############################################################################################

        nijk_remain_df = FilterMinReplicateNumberRunner(
            nijk_remain_df).get_variant_read_count_delete_df(min_replicate_number)
        nijk_remain_df = nijk_remain_df.loc[
            (nijk_remain_df.filter_delete == 0)]
        nijk_remain_df.drop('filter_delete', axis=1, inplace=True)

        # Delete object

        return nijk_remain_df[
            ['run_id', 'marker_id', 'biosample_id', 'variant_id']].copy().drop_duplicates(
            inplace=False)

