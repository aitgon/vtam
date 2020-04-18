# -*- coding: utf-8 -*-
"""LFN Filters

This module will store and run these LFNFilters:

- LFN_var, 2, f2_f4_lfn_delete_variant
- LFN var replicate series, 3, f3_f5_lfn_delete_variant_replicate
- LFN_var_dep , 4, f2_f4_lfn_delete_delete_variant
- LFN vardep_replicate series, 5, f3_f5_lfn_delete_variant_replicate
- LFN_repl, 6, f6_lfn_delete_biosample_replicate_delete
- LFN_readcount, 7, mark_delete_lfn_absolute_read_count
- LFN_all, 8, mark_delete_lfn_do_not_pass_all_filters

Expected results and descriptions are given in the docstrings and in this file:
vtam/discussion_reda_aitor/example_filter.ods

"""
import sys

from vtam.utils.VTAMexception import VTAMexception
from vtam.utils.Logger import Logger
from vtam.utils.VariantReadCountDF import VariantReadCountDF

import pandas


class FilterLFNrunner:

    def __init__(self, variant_read_count_df):
        self.variant_read_count_df = variant_read_count_df[
            ['marker_id', 'run_id', 'variant_id', 'biosample_id', 'replicate', 'read_count']]
        # Instance for all the calculations N_i, N_ij, etc
        self.variant_read_count_lfn_df = VariantReadCountDF(variant_read_count_df)
        #
        if self.variant_read_count_df.shape[1] != 6:
            raise Exception('Columns missing in the variant2sample2replicate2count data frame!')
        #
        ################################
        #  Output variant_read_count_input_df with deleted variants
        ################################
        self.variant_read_count_filter_delete_df = pandas.DataFrame(data={'run_id': [], 'marker_id': [], 'biosample_id': [], 'variant_id': [],
                                                                          'replicate': [],
                                                                          'read_count': [], 'filter_id': [],
                                                                          'filter_delete': []},
                                                                    dtype='uint32')
    def mark_delete_lfn_per_Ni_or_Nik_or_Njk(self, lfn_denominator, threshold, threshold_specific_df=None):
        """

        :param lfn_denominator: string that takes values either: 'N_i', 'N_ik' or 'N_jk'
        :param threshold: float with general threshold
        :param threshold_specific_df: DataFrame with either variant-specific (N_i) or variant-replicate-specific
        deletion threshold
        :return: None: The output of this filter is added to the 'self.variant_read_count_filter_delete_df'
            with filter_id=2 and 'filter_delete'=1 or 0 (General threshold)
            and with filter_id=4 and 'filter_delete'=1 or 0 (Variant-specific threshold)
        """
        if lfn_denominator == 'N_i':  # variant
            this_filter_id = 2
            N_df = self.variant_read_count_lfn_df.get_N_i_df() # Compute N_i_df
            filter_df = self.variant_read_count_df.merge(N_df, left_on=['run_id', 'marker_id', 'variant_id'],
                                                   right_on=['run_id', 'marker_id', 'variant_id'])
            filter_df['lfn_ratio'] = filter_df.read_count / filter_df.N_i
            filter_df['filter_id'] = this_filter_id
            filter_df['threshold'] = threshold
            if not threshold_specific_df is None:
                this_filter_id = 4
                for rowtuple in threshold_specific_df.itertuples():
                    variant_id = rowtuple.variant_id
                    threshold = rowtuple.threshold
                    filter_df['filter_id'] = this_filter_id
                    filter_df.loc[(filter_df.variant_id == variant_id), 'threshold'] = threshold
        elif lfn_denominator == 'N_ik':  # variant_replicate
            this_filter_id = 3
            N_df = self.variant_read_count_lfn_df.get_N_ik_df() # Compute N_ik_df
            filter_df = self.variant_read_count_df.merge(N_df, left_on=['run_id', 'marker_id', 'variant_id', 'replicate'],
                                                   right_on=['run_id', 'marker_id', 'variant_id', 'replicate'])
            filter_df['lfn_ratio'] = filter_df.read_count / filter_df.N_ik
            filter_df['filter_id'] = this_filter_id
            filter_df['threshold'] = threshold
            if not threshold_specific_df is None:
                this_filter_id = 5
                for rowtuple in threshold_specific_df.itertuples():
                    variant_id = rowtuple.variant_id
                    replicate = rowtuple.replicate
                    threshold = rowtuple.threshold
                    filter_df['filter_id'] = this_filter_id
                    filter_df.loc[(filter_df.variant_id == variant_id) & (filter_df.replicate == replicate), 'threshold'] = threshold
        elif lfn_denominator == 'N_jk':  # biosample_replicate
            this_filter_id = 6
            N_df = self.variant_read_count_lfn_df.get_N_jk_df() # Compute N_jk_df
            filter_df = self.variant_read_count_df.merge(N_df, left_on=['run_id', 'marker_id', 'biosample_id', 'replicate'],
                                                   right_on=['run_id', 'marker_id', 'biosample_id', 'replicate'])
            filter_df['lfn_ratio'] = filter_df.read_count / filter_df.N_jk
            filter_df['filter_id'] = this_filter_id
            filter_df['threshold'] = threshold
        else:
            Logger.instance().critical(VTAMexception("Internal error. VTAM will exit."))
            sys.exit(0)
        #
        # Initialize filter: Keep everything
        filter_df['filter_delete'] = False
        #
        # Mark for deletion all variants with read_count=0
        filter_df.loc[
            filter_df.read_count == 0, 'filter_delete'] = True
        #
        # Mark for deletion all filters with 'lfn_ratio'<lfn_variant_threshold
        filter_df.loc[
            filter_df['lfn_ratio'] < filter_df['threshold'], 'filter_delete'] = True

        #
        #  Keep important columns
        filter_df = filter_df[['run_id', 'marker_id', 'biosample_id', 'replicate', 'variant_id', 'read_count',
                   'filter_id', 'filter_delete']]
        #
        # Prepare output variant_read_count_input_df and concatenate vertically output variant_read_count_input_df to self.variant_read_count_filter_delete_df
        self.variant_read_count_filter_delete_df = pandas.concat([self.variant_read_count_filter_delete_df, filter_df],
                                                                 sort=False, axis=0)

    def mark_delete_lfn_absolute_read_count(self, lfn_read_count_threshold):
        """
        Low frequency noise filter (LFN_readcount) with a single threshold lfn_read_count_threshold.
        Function calculating the Low Frequency Noise per users defined minimal readcount
        Function IDs: 7 (LFN_readcount)

        This function implements filter mark_delete_lfn_absolute_read_count (LFN_readcount)
        lfn_read_count_threshold = 10


        This filters deletes the variant if the read count N_ijk of variant i in biosample j
        and replicate k is below threshold lfn_read_count_threshold.
        The deletion condition is: N_ijk < lfn_read_count_threshold .


        Pseudo-algorithm of this function:

           1. Set variant/biosample/replicate for deletion if read_count N_ijk = 0
           2. Set variant/biosample/replicate for deletion if N_ijk < lfn_read_count_threshold



        Updated:
           February 24, 2019

        Args:
           lfn_read_count_threshold (float): Default deletion threshold


        Returns:
           None: The output of this filter is added to the 'self.variant_read_count_filter_delete_df'
           with filter_id='mark_delete_lfn_absolute_read_count' and 'filter_delete'= 1 or 0

        """
        this_filter_id = 7
        # Selecting all the indexes where the ratio is below the minimal readcount
        filter_df = self.variant_read_count_df.copy()

        # do_not_pass_variant_id_list = self.variant_read_count_input_df.loc[
        #     self.variant_read_count_input_df['read_count'] < lfn_read_count_threshold].variant_id.tolist()
        # do_not_pass_replicate_list = self.variant_read_count_input_df.loc[
        #     self.variant_read_count_input_df['read_count'] < lfn_read_count_threshold].replicate.tolist()

        filter_df['filter_id'] = this_filter_id  #  set this filter
        filter_df['filter_delete'] = False

        filter_df.loc[
            filter_df.read_count < lfn_read_count_threshold, 'filter_delete'] = True
        filter_df = filter_df[['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate', 'read_count',
                   'filter_id', 'filter_delete']]
        #
        # Concatenate vertically output variant_read_count_input_df
        #  Prepare output variant_read_count_input_df and concatenate to self.variant_read_count_filter_delete_df

        self.variant_read_count_filter_delete_df = pandas.concat([self.variant_read_count_filter_delete_df, filter_df],
                                                                 sort=False)

    def mark_delete_lfn_do_not_pass_all_filters(self):
        this_filter_id = 8
        # Calculating the total of reads by variant
        filter_df = self.variant_read_count_filter_delete_df[
            ['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate',
             'filter_delete']].groupby(['run_id', 'marker_id', 'variant_id', 'biosample_id',
                                        'replicate']).agg({'filter_delete': sum})
        filter_df = filter_df.rename(columns={'filter_delete': 'filter_delete_aggregate'})

        # Merge the column with the total reads by sample replicates for calculate the ratio
        filter_df = self.variant_read_count_df.merge(filter_df, left_on=['run_id', 'marker_id', 'variant_id', 'biosample_id',
                                                             'replicate'],
                                               right_on=['run_id', 'marker_id', 'variant_id', 'biosample_id',
                                                         'replicate'])
        #
        # Initialize
        filter_df['filter_id'] = this_filter_id
        filter_df['filter_delete'] = False
        #
        #
        # Selecting all the indexes where the ratio is below the ratio
        filter_df.loc[
            filter_df.filter_delete_aggregate > 0, 'filter_delete'] = True
        #
        # Let only interest columns
        filter_df = filter_df[['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate', 'read_count',
                   'filter_id', 'filter_delete']]
        #
        self.variant_read_count_filter_delete_df = pandas.concat([self.variant_read_count_filter_delete_df, filter_df],
                                                                 sort=False)
