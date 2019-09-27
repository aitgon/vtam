# -*- coding: utf-8 -*-
"""LFN Filters

This module will store and run these LFNFilters:

- LFN_var, 2, f2_f4_lfn_delete_variant
- LFN var replicate series, 3, f3_f5_lfn_delete_variant_replicate
- LFN_var_dep , 4, f2_f4_lfn_delete_delete_variant
- LFN vardep_replicate series, 5, f3_f5_lfn_delete_variant_replicate
- LFN_repl, 6, f6_lfn_delete_biosample_replicate_delete
- LFN_readcount, 7, f7_lfn_delete_absolute_read_count
- LFN_all, 8, f8_lfn_delete_do_not_pass_all_filters

Expected results and descriptions are given in the docstrings and in this file:
vtam/discussion_reda_aitor/example_filter.ods

"""

from wopmetabarcoding.utils.logger import logger
import inspect
import pandas


def f1_lfn_delete_singleton(run_biosample_replicate_variant_read_count_df):
    """
    This function removes variants that appear just once over all sample/replicates


    Args:
        run_biosample_replicate_variant_read_count_df (Pandas Dataframe): Dataframe with this columns: run_id, variant_sequence, biosample_id, replicate_id, read_count

    Returns: Pandas DataFrame with same columns as input, but without singletons (Variant with just one read)

    """
    df = run_biosample_replicate_variant_read_count_df[['variant_sequence', 'read_count']].groupby(
        by=['variant_sequence']).sum().reset_index()

    # Get names of indexes for which variant_sequence has the read_count equal 1
    singleton_sequence_list = df[df['read_count'] == 1].variant_sequence.tolist()

    # Delete these row indexes from dataFrame
    passed_lfn_delete_singleton_df = run_biosample_replicate_variant_read_count_df.loc[
        ~run_biosample_replicate_variant_read_count_df.variant_sequence.isin(singleton_sequence_list),]
    return passed_lfn_delete_singleton_df



class FilterLFNRunner:


    def __init__(self, variant_read_count_df):
        logger.debug(
            "file: {}; line: {}; FilterNonLFNRunner.__init__".format(__file__, inspect.currentframe().f_lineno))
        # self.variant_df = variant_df
        self.variant_read_count_df = variant_read_count_df[['marker_id', 'run_id', 'variant_id', 'biosample_id', 'replicate_id', 'read_count']]
        # self.marker_id = marker_id
        #
        if self.variant_read_count_df.shape[1] != 6:
            raise Exception('Columns missing in the variant2sample2replicate2count data frame!')
        #
        ################################
        # Output df with deleted variants
        ################################
        self.delete_variant_df = pandas.DataFrame(data={'variant_id':[], 'biosample_id':[], 'replicate_id':[], 'read_count':[], 'filter_id':[], 'filter_delete':[]}, dtype='int')
        logger.debug(
            "file: {}; line: {}; Initial nb of variants {}".format(__file__, inspect.currentframe().f_lineno,
                                                               (self.delete_variant_df.sum(axis=1) == self.delete_variant_df.shape[1]).sum()))


    def f2_f4_lfn_delete_variant(self, lfn_variant_threshold, threshold_specific_df=None):
        """
        This filter deletes the variant i in biosample j and replicate k
        if the ratio of its read count N_ijk to the total read_count N_i of variant i
        is below a threshold parameter lfn_per_variant_threshold.
        Function IDs: 2 (lfn_var_threshold_specific is None) ou 4 (lfn_var_threshold_specific is not None)

        The argument lfn_var_threshold_specific allows a dictionary like {9: 0.05, 22: 0.01} with a variant-specific
        threshold.

        Pseudo-algorithm of this function:

        1. Compute ratio N_ijk / N_i
        2. Set variant/biosample/replicate for deletion if read_count N_ijk = 0
        3. Set variant/biosample/replicate for deletion if ratio N_ijk / N_i < lfn_per_variant_threshold
        4. If variant specific thresholds, copy these variant/biosample/replicate rows
          4.1 Set variant/biosample/replicate for deletion if N_ijk / N_i < lfn_var_threshold_specific_i

        Updated:
        May 28, 2019

        Args:
            lfn_variant_threshold (float): Default deletion threshold
            threshold_specific_df (:obj:`dict`, optional): Variant-specific deletion threshold

        Returns:
            None: The output of this filter is added to the 'self.delete_variant_df'
            with filter_id=2 and 'filter_delete'=1 or 0 (General threshold)
            and with filter_id=4 and 'filter_delete'=1 or 0 (Variant-specific threshold)


        """
        this_filter_id = 2
        # Write log
        logger.debug(
            "file: {}; line: {}; {}".format(__file__, inspect.currentframe().f_lineno, this_filter_id))
        ######################
        # Calculating the total of reads by variant
        df2 = self.variant_read_count_df.copy()
        df2 = df2[['run_id', 'marker_id', 'variant_id', 'read_count']].groupby(by=['run_id', 'marker_id', 'variant_id']).sum().reset_index()
        # Merge the column with the total reads by variant for calculate the ratio
        df2 = self.variant_read_count_df.merge(df2, left_on=['run_id', 'marker_id', 'variant_id'], right_on=['run_id', 'marker_id', 'variant_id'])
        df2 = df2.rename(columns={'read_count_x': 'read_count'})
        df2 = df2.rename(columns={'read_count_y': 'read_count_per_variant'})
        # Calculate the ratio
        df2['low_frequence_noice_per_variant'] = df2.read_count / df2.read_count_per_variant
        #
        # Initialize filter: Keep everything
        df2['filter_id'] = this_filter_id
        df2['filter_delete'] = False
        #
        # Mark for deletion all filters with read_count=0
        df2.loc[
            df2.read_count == 0, 'filter_delete'] = True
        #
        # Mark for deletion all filters with low_frequence_noice_per_variant<lfn_variant_threshold
        df2.loc[
            df2.low_frequence_noice_per_variant < lfn_variant_threshold, 'filter_delete'] = True

        ####################################################
        # LFN_var_dep if there are variant specific thresholds: threshold_specific_df not None
        ####################################################
        if not threshold_specific_df is None:
            this_filter_id = 4
            for rowtuple in threshold_specific_df.itertuples():
                variant_id = rowtuple.variant_id
                threshold = rowtuple.threshold

                #
                df2_f4_variant_id = df2.loc[(df2.variant_id == variant_id)].copy()
                #
                # Initialize filter: Keep everything
                df2_f4_variant_id.loc[df2_f4_variant_id.variant_id == variant_id, 'filter_id'] = this_filter_id
                df2_f4_variant_id.loc[df2_f4_variant_id.variant_id == variant_id, 'filter_delete'] = False
                #
                # Mark for deletion all filters with low_frequence_noice_per_variant<lfn_variant_threshold
                df2_f4_variant_id.loc[
                    df2_f4_variant_id.low_frequence_noice_per_variant < threshold, 'filter_delete'] = True

                #
                # Keep important columns
                df2_f4_variant_id = df2_f4_variant_id[['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id', 'read_count',
                           'filter_id', 'filter_delete']]
                #
                # Concatenate vertically output df
                # Prepare output df and concatenate to self.delete_variant_df
                self.delete_variant_df = pandas.concat([self.delete_variant_df, df2_f4_variant_id], sort=False)

        #
        # Keep important columns
        df2 = df2[['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id', 'read_count',
                   'filter_id', 'filter_delete']]
        # Concatenate vertically output df
        # Prepare output df and concatenate to self.delete_variant_df
        self.delete_variant_df = pandas.concat([self.delete_variant_df, df2], sort=False)


    def f3_f5_lfn_delete_variant_replicate(self, lfn_variant_replicate_threshold, threshold_specific_df=None):
        """
        This filter deletes the variant i in biosample j and replicate k
        if the ratio of its read count N_ijk to the replicate read_count N_ik
        is below a threshold parameter lfn_per_replicate_series_threshold.
        Function IDs: 3 (lfn_per_replicate_threshold_specific is None) ou 5 (lfn_per_replicate_threshold_specific is not None).

        Pseudo-algorithm of this function:

            1. Compute ratio N_ijk / N_ik
            2. Set variant/biosample/replicate for deletion if read_count N_ijk = 0
            3. Set variant/biosample/replicate for deletion if ratio N_ijk / N_ik < lfn_per_replicate_series_threshold
            4. If variant specific thresholds, test the ratio of these variant/biosample/replicate rows
             4.1 Set variant/biosample/replicate for deletion if N_ijk / N_ik < lfn_var_threshold_specific_i

         Updated:
            May 28, 2019

        Args:
            lfn_variant_replicate_threshold (float): Default deletion threshold
            threshold_specific_df (:obj:`dict`, optional): Variant-specific deletion threshold


        Returns:
            None: The output of this filter is added to the 'self.delete_variant_df'
            with filter_id=3 and 'filter_delete'=1 or 0 (General threshold)
            or with filter_id=5 and 'filter_delete'=1 or 0 (Variant-specific threshold)

        """
        this_filter_id = 3
        # Get function
        # Write log
        logger.debug(
            "file: {}; line: {}; {}".format(__file__, inspect.currentframe().f_lineno, this_filter_id))
        ######################
        df2 = self.variant_read_count_df.copy()
        # Calculating the total of reads by replicate series
        df2 = df2[['run_id', 'marker_id', 'variant_id', 'replicate_id', 'read_count']].groupby(by=['run_id', 'marker_id', 'variant_id', 'replicate_id']).sum().reset_index()
        # Merge the column with the total reads by variant for calculate the ratio

        df2 = self.variant_read_count_df.merge(df2, left_on=['run_id', 'marker_id', 'variant_id', 'replicate_id'], right_on=['run_id', 'marker_id', 'variant_id', 'replicate_id'])
        df2 = df2.rename(columns={'read_count_x': 'read_count'})
        df2 = df2.rename(columns={'read_count_y': 'read_count_per_replicate_series'})

        # Calculate the ratio
        df2['low_frequence_noice_per_replicate_series'] = df2.read_count / df2.read_count_per_replicate_series

        #
        df2['filter_id'] = this_filter_id # set this filter
        df2['filter_delete'] = False
        df2.loc[
            df2.low_frequence_noice_per_replicate_series < lfn_variant_replicate_threshold, 'filter_delete'] = True
        if not threshold_specific_df is None:
            this_filter_id = 5
            # for variant_id in threshold_specific_df:
            for rowtuple in threshold_specific_df.itertuples():
                variant_id = rowtuple.variant_id
                replicate_id = rowtuple.replicate_id
                threshold = rowtuple.threshold
                # threshold = threshold_specific_df[variant_id]
                df2_f6_variant_id = df2.loc[((df2.variant_id == variant_id) & (df2.replicate_id == replicate_id))].copy()

                # Initialize filter: Keep everything
                df2_f6_variant_id.loc[df2_f6_variant_id.variant_id == variant_id, 'filter_id'] = this_filter_id
                df2_f6_variant_id.loc[df2_f6_variant_id.variant_id == variant_id, 'filter_delete'] = False

                # Mark for deletion all filters with low_frequence_noice_per_variant<lfn_per_variant_threshold
                df2_f6_variant_id.loc[
                    df2_f6_variant_id.low_frequence_noice_per_replicate_series < threshold, 'filter_delete'] = True

                #  Keep important columns
                df2_f6_variant_id = df2_f6_variant_id[['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id', 'read_count',
                                                       'filter_id', 'filter_delete']]

                # Concatenate vertically output df
                #  Prepare output df and concatenate to self.delete_variant_df
                self.delete_variant_df = pandas.concat([self.delete_variant_df, df2_f6_variant_id], sort=False)
        df2 = df2[['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id', 'read_count',
                   'filter_id', 'filter_delete']]

      # Concatenate vertically output df
      # Prepare output df and concatenate to self.delete_variant_df
        self.delete_variant_df = pandas.concat([self.delete_variant_df, df2], sort=False)

    def f6_lfn_delete_biosample_replicate(self, lfn_biosample_replicate_threshold):
        """
        This filter deletes the variant i in biosample j and replicate k
        if the ratio of its read count N_ijk to the read_count N_jk
        is below a threshold parameter lfn_biosample_replicate_threshold (Default 0.001).
        Function IDs: 6 (LFN_repl)

        Pseudo-algorithm of this function:

        1. Compute ratio N_ijk / N_jk
        2. Set variant/biosample/replicate for deletion if read_count N_ijk = 0
        3. Set variant/biosample/replicate for deletion if ratio N_ijk / N_jk < lfn_per_biosample_per_replicate_threshold

        Updated:
            February 23, 2019

        Args:
            lfn_biosample_replicate_threshold (float): Default deletion threshold


        Returns:
            None: The output of this filter is added to the 'self.delete_variant_df'
            with filter_id='f6_lfn_delete_biosample_replicate'= 1 or 0
        """

        this_filter_id = 6
        logger.debug(
            "file: {}; line: {}; {}".format(__file__, inspect.currentframe().f_lineno, this_filter_id))
        # Calculating the total number of reads by sample replicates

        df2 = self.variant_read_count_df[['run_id', 'marker_id', 'biosample_id', 'replicate_id', 'read_count']].groupby(
            by=['run_id', 'marker_id', 'biosample_id', 'replicate_id']).sum().reset_index()


        # Merge the column with the total reads by sample replicates for calculate the ratio
        df2 = self.variant_read_count_df.merge(df2, left_on=['run_id', 'marker_id', 'biosample_id', 'replicate_id'],
                                               right_on=['run_id', 'marker_id', 'biosample_id', 'replicate_id'])
        df2 = df2.rename(columns={'read_count_x': 'read_count'})
        df2 = df2.rename(columns={'read_count_y': 'read_count_per_biosample_replicate'})
        #
        # Initialize
        this_filter_name = 'f6_lfn_delete_biosample_replicate'
        df2['filter_id'] = this_filter_id
        df2['filter_delete'] = False
        #
        # Calculate the ratio
        df2[
            'low_frequence_noice_per_replicate'] = df2.read_count / df2.read_count_per_biosample_replicate
        #
        # Selecting all the indexes where the ratio is below the ratio
        df2.loc[
            df2.low_frequence_noice_per_replicate < lfn_biosample_replicate_threshold , 'filter_delete'] = True
        #
        #Let only interest columns
        df2 = df2[['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id', 'read_count',
                       'filter_id', 'filter_delete']]
        #
        self.delete_variant_df = pandas.concat([self.delete_variant_df, df2], sort=False)


    def f7_lfn_delete_absolute_read_count(self, lfn_read_count_threshold):
        """
        Low frequency noise filter (LFN_readcount) with a single threshold lfn_read_count_threshold.
        Function calculating the Low Frequency Noise per users defined minimal readcount
        Function IDs: 7 (LFN_readcount)

        This function implements filter f7_lfn_delete_absolute_read_count (LFN_readcount)
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
           None: The output of this filter is added to the 'self.delete_variant_df'
           with filter_id='f7_lfn_delete_absolute_read_count' and 'filter_delete'= 1 or 0

        """
        this_filter_id = 7
        logger.debug(
            "file: {}; line: {}; {}".format(__file__, inspect.currentframe().f_lineno, this_filter_id))
        # Selecting all the indexes where the ratio is below the minimal readcount
        df2 = self.variant_read_count_df.copy()
        do_not_pass_variant_id_list = self.variant_read_count_df.loc[
            self.variant_read_count_df['read_count'] < lfn_read_count_threshold].variant_id.tolist()
        do_not_pass_replicate_id_list = self.variant_read_count_df.loc[
            self.variant_read_count_df['read_count'] < lfn_read_count_threshold].replicate_id.tolist()

        df2['filter_id'] = this_filter_id  #  set this filter
        df2['filter_delete'] = False

        df2.loc[
            df2.read_count < lfn_read_count_threshold, 'filter_delete'] = True
        df2 = df2[['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id', 'read_count',
                   'filter_id', 'filter_delete']]
        #
        # Concatenate vertically output df
        #  Prepare output df and concatenate to self.delete_variant_df


        self.delete_variant_df = pandas.concat([self.delete_variant_df, df2], sort=False)
        logger.debug(
            "file: {}; line: {}; Nb variants passed {}".format(__file__, inspect.currentframe().f_lineno,
                                                               (self.delete_variant_df.sum(axis=1) ==
                                                                self.delete_variant_df.shape[1]).sum()))



    def f8_lfn_delete_do_not_pass_all_filters(self):
        this_filter_id = 8
        # Calculating the total of reads by variant
        df2 = self.delete_variant_df[['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id',
                                      'filter_delete']].groupby(['run_id', 'marker_id', 'variant_id', 'biosample_id',
                                                                 'replicate_id']).sum().reset_index()
        df2 = df2.rename(columns={'filter_delete': 'filter_delete_aggregate'})

        # Merge the column with the total reads by sample replicates for calculate the ratio
        df2 = self.variant_read_count_df.merge(df2, left_on=['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id'],
                                               right_on=['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id'])
        #
        # Initialize
        df2['filter_id'] = this_filter_id
        df2['filter_delete'] = False
        #
        #
        # Selecting all the indexes where the ratio is below the ratio
        df2.loc[
            df2.filter_delete_aggregate>0, 'filter_delete'] = True
        #
        #Let only interest columns
        df2 = df2[['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id', 'read_count',
                       'filter_id', 'filter_delete']]
        #
        self.delete_variant_df = pandas.concat([self.delete_variant_df, df2], sort=False)

