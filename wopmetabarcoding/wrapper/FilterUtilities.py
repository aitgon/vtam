import inspect

from sqlalchemy import select
from wopmetabarcoding.utils.PathFinder import PathFinder
from wopmetabarcoding.utils.VSearch import VSearch1, Vsearch2, Vsearch3
import pandas, itertools
from Bio import SeqIO
import re
import os
from wopmetabarcoding.utils.constants import tempdir

from math import floor

from wopmetabarcoding.utils.logger import logger


class FilterRunner:

    def __init__(self, variant_df, variant_read_count_df, marker_id):
        logger.debug(
            "file: {}; line: {}; FilterRunner.__init__".format(__file__, inspect.currentframe().f_lineno))
        self.variant_df = variant_df
        self.variant_read_count_df = variant_read_count_df
        self.marker_id = marker_id
        #
        self.tempdir = os.path.join(tempdir, "FilterUtilities", "marker_id_{}".format(marker_id), "FilterUtilities", self.__class__.__name__)
        PathFinder.mkdir_p(self.tempdir)
        #
        if self.variant_read_count_df.shape[1] != 4:
            raise Exception('Columns missing in the variant2sample2replicate2count data frame!')
        #
        # Empty table: first time of executing this code
        #
        # passed_variant_df
        filters = [
            'passed',
            'f1_lfn1_per_replicate',
            'f2_lfn2_per_variant',
            'f3_lfn2_per_replicate_series',
            'f4_lfn3_read_count',
            'f5_lfn4_per_variant_with_cutoff',
            'f6_lfn4_per_replicate_series_with_cutoff',
            'f7_min_repln',
            'f8_min_replp',
            'f9_pcr_error',
            'f10_chimera',
            'f10_chimera_borderline',
            'f11_renkonen',
            'f12_indel',
        ]
        #
        ############################################
        # DF with results of filters. This DF is composed of
        #
        # Col 1: 'variant_id' (integer, index, unique and primary key)
        # Col 2: 'passed'. Did this variant_id passed all the filters (Boolean, default True)
        # Col 3: 'f1_lfn1_per_replicate'. Did this variant_id passed filter 'f1_lfn1_per_replicate' (Boolean, default True)
        # Col 4: ...
        # Col 15: 'indel'. Did this variant_id passed filter 'indel' (Boolean, default True)
        ############################################
        # 
        # Create
        # Default state: variant passed filters
        self.passed_variant_df = pandas.DataFrame(dict(zip(filters, itertools.repeat(True))), index=self.variant_df.id.unique().tolist())
        ################################
        # Fix output avec Reda Mekdad:
        ################################
        self.passed_variant_mekdad_df = pandas.DataFrame(data={'variant_id':[], 'biosample_id':[], 'replicate_id':[], 'filter_name':[], 'filter_passed':[]})
        logger.debug(
            "file: {}; line: {}; Initial nb of variants {}".format(__file__, inspect.currentframe().f_lineno,
                                                               (self.passed_variant_df.sum(axis=1) == self.passed_variant_df.shape[1]).sum()))

    def f1_lfn1_per_replicate(self, lfn_per_replicate_threshold):
        """
        Function returning indices passing this filter.
        It calculates the Low Frequency Noise per sample replicate
        :param variant_read_count_df: dataframe containing the information
        :param lfn_per_replicate_threshold: threshold defined by the user
        :return: List of the index which don't pass the filter
        """
        this_filter_name = inspect.stack()[0][3]
        logger.debug(
            "file: {}; line: {}; {}".format(__file__, inspect.currentframe().f_lineno, this_filter_name))
        # Calculating the total number of reads by sample replicates
        df2 = self.variant_read_count_df[['biosample_id', 'replicate_id', 'read_count']].groupby(by=['biosample_id', 'replicate_id']).sum().reset_index()
        # Merge the column with the total reads by sample replicates for calculate the ratio
        df2 = self.variant_read_count_df.merge(df2, left_on=['biosample_id', 'replicate_id'], right_on=['biosample_id', 'replicate_id'])
        df2.columns = ['variant_id', 'biosample_id', 'replicate_id', 'read_count_per_variant_per_biosample_per_replicate', 'read_count_per_biosample_replicate']
        # Calculate the ratio
        df2['low_frequence_noice_per_replicate'] = df2.read_count_per_variant_per_biosample_per_replicate / df2.read_count_per_biosample_replicate
        # Selecting all the indexes where the ratio is below the ratio
        do_not_pass_variant_id_list = df2.loc[df2.low_frequence_noice_per_replicate < lfn_per_replicate_threshold,'variant_id'].tolist()
        #
        self.passed_variant_df.loc[do_not_pass_variant_id_list, 'passed'] = False
        self.passed_variant_df.loc[do_not_pass_variant_id_list, this_filter_name] = False
        logger.debug(
            "file: {}; line: {}; Nb variants passed {}".format(__file__, inspect.currentframe().f_lineno,
                                                               (self.passed_variant_df.sum(axis=1) == self.passed_variant_df.shape[1]).sum()))


    def f2_lfn2_per_variant(self, lfn_per_variant_threshold):
        """
        This filter corresponds to LFN_var example of Emese

        :param variant_read_count_df: dataframe containing the information
        :param lfn_per_variant_threshold: threshold defined by the user
        :return: List of the index which don't pass the filter
        """
        this_filter_name = inspect.stack()[0][3]
        logger.debug(
            "file: {}; line: {}; {}".format(__file__, inspect.currentframe().f_lineno, this_filter_name))
        #########################
        # Calculating the total of reads by variant
        df2 = self.variant_read_count_df[['variant_id', 'read_count']].groupby(by=['variant_id']).sum().reset_index()
        # Merge the column with the total reads by variant for calculate the ratio
        df2 = self.variant_read_count_df.merge(df2, left_on='variant_id', right_on='variant_id')
        df2.columns = ['variant_id', 'biosample_id', 'replicate_id', 'read_count_per_variant_per_biosample_per_replicate', 'read_count_per_variant']
        # Calculate the ratio
        df2['low_frequence_noice_per_variant'] = df2.read_count_per_variant_per_biosample_per_replicate / df2.read_count_per_variant
        # Selecting all the indexes where the ratio is below the ratio
        do_not_pass_variant_id_list = df2.loc[df2.low_frequence_noice_per_variant < lfn_per_variant_threshold].variant_id.tolist() # select rare reads to discard later
        #
        self.passed_variant_df.loc[do_not_pass_variant_id_list, 'passed'] = False
        self.passed_variant_df.loc[do_not_pass_variant_id_list, this_filter_name] = False
        # import pdb; pdb.set_trace()
        logger.debug(
            "file: {}; line: {}; Nb variants passed {}".format(__file__, inspect.currentframe().f_lineno,
                                                               (self.passed_variant_df.sum(axis=1) == self.passed_variant_df.shape[1]).sum()))


    def f2_lfn2_per_variant_mekdad(self, lfn_per_variant_threshold, lfn_per_variant_threshold_i=None):
        """
        This filter corresponds to LFN_var example of Emese

        :param variant_read_count_df: dataframe containing the information
        :param lfn_per_variant_threshold: threshold defined by the user
        :return: None. Result is in 'f2_lfn2_per_variant_mekda' column. True if variant-biosample-replicated passed
        the filter or False otherwise
        """
        this_filter_name = inspect.stack()[0][3] # Get function
        # Write log
        logger.debug(
            "file: {}; line: {}; {}".format(__file__, inspect.currentframe().f_lineno, this_filter_name))
        ######################
        # Calculating the total of reads by variant
        df2 = self.variant_read_count_df[['variant_id', 'read_count']].groupby(by=['variant_id']).sum().reset_index()
        # Merge the column with the total reads by variant for calculate the ratio
        df2 = self.variant_read_count_df.merge(df2, left_on='variant_id', right_on='variant_id')
        df2 = df2.rename(columns={'read_count_x': 'read_count_per_variant_per_biosample_per_replicate'})
        df2 = df2.rename(columns={'read_count_y': 'read_count_per_variant'})
        # Calculate the ratio
        df2['low_frequence_noice_per_variant'] = df2.read_count_per_variant_per_biosample_per_replicate / df2.read_count_per_variant
        #
        df2['filter_name'] = this_filter_name
        df2['filter_passed'] = True
        # df2[this_filter_name] = True
        df2.loc[
            df2.low_frequence_noice_per_variant < lfn_per_variant_threshold, 'filter_passed'] = False
        df2 = df2[['variant_id', 'biosample_id', 'replicate_id',
                   'filter_name', 'filter_passed']]
        #
        # Concatenate vertically output df
        # Prepare output df and concatenate to self.passed_variant_mekdad_df
        self.passed_variant_mekdad_df = pandas.concat([self.passed_variant_mekdad_df, df2], sort=False)


    def f3_lfn2_per_replicate_series_mekdad(self, lfn_per_replicate_series_threshold, lfn_per_variant_threshold_i=None):
        """
        This filter corresponds to LFN_var example of Emese

        :param variant_read_count_df: dataframe containing the information
        :param lfn_per_variant_threshold: threshold defined by the user
        :return: None. Result is in 'f2_lfn2_per_variant_mekda' column. True if variant-biosample-replicated passed
        the filter or False otherwise
        """
        this_filter_name = inspect.stack()[0][3] # Get function
        # Write log
        logger.debug(
            "file: {}; line: {}; {}".format(__file__, inspect.currentframe().f_lineno, this_filter_name))
        ######################
        # Calculating the total of reads by variant
        df2 = self.variant_read_count_df[['replicate_id', 'read_count']].groupby(by=['replicate_id']).sum().reset_index()
        # Merge the column with the total reads by replicate for calculate the ratio
        df2 = self.variant_read_count_df.merge(df2, left_on='replicate_id', right_on='replicate_id')
        df2 = df2.rename(columns={'read_count_x': 'read_count_per_variant_per_biosample_replicate'})
        df2 = df2.rename(columns={'read_count_y': 'read_count_per_replicate_series'})
        # Calculate the ratio
        df2['low_frequence_noice_per_replicate_series'] = df2.read_count_per_variant_per_biosample_replicate / df2.read_count_per_replicate_series
        # ####################
        # Output
        self.passed_variant_mekdad_df[this_filter_name] = True
        self.passed_variant_mekdad_df.loc[
            df2.low_frequence_noice_per_replicate_series < lfn_per_replicate_series_threshold, this_filter_name] = False

    def f3_lfn2_per_replicate_series(self, lfn_per_replicate_series_threshold):
        """
        Function calculating the Low Frequency Noise per replicate series
        :param variant_read_count_df: dataframe containing the information
        :param lfn_per_replicate_series_threshold: threshold defined by the user
        :return: List of the index which don't pass the filter
        """
        this_filter_name = inspect.stack()[0][3]
        logger.debug(
            "file: {}; line: {}; {}".format(__file__, inspect.currentframe().f_lineno, this_filter_name))
        # Calculating the total of reads by replicate series
        df2 = self.variant_read_count_df[['replicate_id', 'read_count']].groupby(by=['replicate_id']).sum().reset_index()
        # Merge the column with the total reads by replicate series for calculate the ratio
        df2 = self.variant_read_count_df.merge(df2, left_on='replicate_id', right_on='replicate_id')
        df2.columns = ['variant_id', 'biosample_id', 'replicate_id',
                       'read_count_per_variant_per_biosample_replicate', 'read_count_per_replicate_series']
        df2[
            'low_frequence_noice_per_replicate_series'] = df2.read_count_per_variant_per_biosample_replicate / df2.read_count_per_replicate_series
        # Selecting all the indexes where the ratio is below the ratio
        do_not_pass_variant_id_list = df2.loc[
                        df2.low_frequence_noice_per_replicate_series < lfn_per_replicate_series_threshold].variant_id.tolist()  #  select rare reads to discard later
        #
        self.passed_variant_df.loc[do_not_pass_variant_id_list, 'passed'] = False
        self.passed_variant_df.loc[do_not_pass_variant_id_list, this_filter_name] = False
        logger.debug(
            "file: {}; line: {}; Nb variants passed {}".format(__file__, inspect.currentframe().f_lineno,
                                                               (self.passed_variant_df.sum(axis=1) == self.passed_variant_df.shape[1]).sum()))

    def f4_lfn3_read_count(self, lfn_read_count_threshold):
        """
        Function calculating the Low Frequency Noise per users defined minimal readcount
        
        :param variant_read_count_df: dataframe containing the information
        :param lfn_read_count_threshold: threshold defined by the user
        :return: List of the index which don't pass the filter
        """
        this_filter_name = inspect.stack()[0][3]
        logger.debug(
            "file: {}; line: {}; {}".format(__file__, inspect.currentframe().f_lineno, this_filter_name))
        # Selecting all the indexes where the ratio is below the minimal readcount
        do_not_pass_variant_id_list = self.variant_read_count_df.loc[self.variant_read_count_df['read_count'] < lfn_read_count_threshold].variant_id.tolist()
        #
        self.passed_variant_df.loc[do_not_pass_variant_id_list, 'passed'] = False
        self.passed_variant_df.loc[do_not_pass_variant_id_list, this_filter_name] = False
        logger.debug(
            "file: {}; line: {}; Nb variants passed {}".format(__file__, inspect.currentframe().f_lineno,
                                                               (self.passed_variant_df.sum(axis=1) == self.passed_variant_df.shape[1]).sum()))

    def f5_lfn4_per_variant_with_cutoff(self, cutoff_tsv):
        """
        Function calculating the Low Frequency Noise per variant against cutoff

        :param variant_read_count_df: dataframe containing the information
        :param cutoff_tsv: file containing the cutoffs for each variant
        :return: List of the index which don't pass the filter
        """
        # TODO: Must be updated for new FilterRunner class
        # TODO: Variant sequence must be taken from variant_df
        raise NotImplementedError("Must be updated for new FilterRunner class")
        this_filter_name = inspect.stack()[0][3]
        logger.debug(
            "file: {}; line: {}; {}".format(__file__, inspect.currentframe().f_lineno, this_filter_name))
        cutoff_df = pandas.read_csv(cutoff_tsv, sep="\t", header=0)
        cutoff_df.columns = ['variant_sequence', 'cutoff']
        cutoff_df['variant_id'] = None
        for variant_sequence in cutoff_df.variant_sequence.tolist():
            with engine.connect() as conn:
                stmt = select([variant_model.__table__.c.id]).where(variant_model.__table__.c.sequence == variant_sequence)
                try:
                    variant_id = conn.execute(stmt).first()[0]
                    cutoff_df.loc[cutoff_df.variant_sequence == variant_sequence, 'variant_id'] = variant_id
                except TypeError:
                    pass
        cutoff_df = cutoff_df.loc[~cutoff_df.variant_id.isna(), ['variant_id', 'cutoff']] # remove absent variants
        #
        # Calculating the total of reads by variant
        df2 = self.variant_read_count_df[['variant_id', 'read_count']].groupby(by=['variant_id']).sum().reset_index()
        # Merge the column with the total reads by variant for calculate the ratio
        df2 = self.variant_read_count_df.merge(df2, left_on='variant_id', right_on='variant_id')
        df2.columns = ['variant_id', 'biosample_id', 'replicate_id', 'read_count_per_variant_per_biosample_replicate', 'read_count_per_variant']
        df2[this_filter_name] = df2.read_count_per_variant_per_biosample_replicate / df2.read_count_per_variant
        #
        # merge with cutoff
        cutoff_df['variant_id'] = pandas.to_numeric(cutoff_df.variant_id)
        df2 = df2.merge(cutoff_df, left_on='variant_id', right_on='variant_id')
        do_not_pass_variant_id_list = df2.ix[df2[this_filter_name] < df2.cutoff].variant_id.tolist()
        #
        self.passed_variant_df.loc[do_not_pass_variant_id_list, 'passed'] = False
        self.passed_variant_df.loc[do_not_pass_variant_id_list, this_filter_name] = False
        logger.debug(
            "file: {}; line: {}; Nb variants passed {}".format(__file__, inspect.currentframe().f_lineno,
                                                               (self.passed_variant_df.sum(axis=1) == self.passed_variant_df.shape[1]).sum()))

    def f6_lfn4_per_replicate_series_with_cutoff(self, cutoff_tsv):
        """
            Function calculating the Low Frequency Noise per replicate series against cutoff

            :param cutoff_tsv: file containing the cutoffs for each variant
            :return: None
            """
        this_filter_name = inspect.stack()[0][3]
        logger.debug(
            "file: {}; line: {}; {}".format(__file__, inspect.currentframe().f_lineno, this_filter_name))
        cutoff_df = pandas.read_csv(cutoff_tsv, sep="\t", header=0)
        cutoff_df.columns = ['variant_sequence', 'cutoff']
        #
        df2 = self.variant_read_count_df.groupby(by=['replicate']).sum()
        df2 = self.variant_read_count_df.merge(df2, left_on='replicate', right_index=True)
        df2.columns = ['variant_sequence', 'replicate', 'biosample', 'sample_replicate',
                       'read_count_per_variant_per_biosample_replicate', 'read_count_per_replicate_series']
        df2[
            'low_frequence_noice_per_replicate_series'] = df2.read_count_per_variant_per_biosample_replicate / df2.read_count_per_replicate_series
        #
        # merge with cutoff
        df2 = df2.merge(cutoff_df, left_on='variant_sequence', right_on='variant_sequence')
        do_not_pass_variant_id_list = df2.ix[df2.low_frequence_noice_per_replicate_series < df2.cutoff].variant_id.tolist()
        #
        self.passed_variant_df.loc[do_not_pass_variant_id_list, 'passed'] = False
        self.passed_variant_df.loc[do_not_pass_variant_id_list, this_filter_name] = False
        logger.debug(
            "file: {}; line: {}; Nb variants passed {}".format(__file__, inspect.currentframe().f_lineno,
                                                               (self.passed_variant_df.sum(axis=1) == self.passed_variant_df.shape[1]).sum()))

    # def execute_filters(self):
    #     indices_to_drop = self.passsed_variant_ids
    #     self.variant_read_count_df.drop(indices_to_drop, inplace=True)
    #     self.passsed_variant_ids = [] # reset passsed_variant_ids

    def f7_min_repln(self, min_count):
        """
        Filter watching if a variant is present the minimal number of sample replicates of a sample

        :param min_count: minimal number of replicates in which the variant must be present
        :return: None
        """
        this_filter_name = inspect.stack()[0][3]
        logger.debug(
            "file: {}; line: {}; {}".format(__file__, inspect.currentframe().f_lineno, this_filter_name))
        for biosample_id in self.variant_read_count_df.biosample_id.unique():
            df_biosample = self.variant_read_count_df.loc[self.variant_read_count_df['biosample_id'] == biosample_id]
            df3 = df_biosample['variant_id'].value_counts().to_frame()
            df3.columns = ['variant_count']
            df_biosample = df_biosample.merge(df3, left_on='variant_id', right_index=True)
            do_not_pass_variant_id_list = df_biosample.loc[
                df_biosample.variant_count < min_count].variant_id.unique().tolist()
            #
            self.passed_variant_df.loc[do_not_pass_variant_id_list, 'passed'] = False
            self.passed_variant_df.loc[do_not_pass_variant_id_list, this_filter_name] = False
        logger.debug(
            "file: {}; line: {}; Nb variants passed {}".format(__file__, inspect.currentframe().f_lineno,
                                                           (self.passed_variant_df.sum(axis=1) == self.passed_variant_df.shape[1]).sum()))


    def f8_min_replp(self):
        """
        Filter watching if a variant is present the more more than 1/3 of the number of sample replicates of a sample

        :return: None
        """
        this_filter_name = inspect.stack()[0][3]
        logger.debug(
            "file: {}; line: {}; {}".format(__file__, inspect.currentframe().f_lineno, this_filter_name))
        for biosample_id in self.variant_read_count_df.biosample_id.unique():
            df_biosample = self.variant_read_count_df.loc[self.variant_read_count_df['biosample_id'] == biosample_id]
            replicate_count = len(df_biosample.replicate_id.unique().tolist())
            df3 = df_biosample['variant_id'].value_counts().to_frame()
            df3.columns = ['variant_count']
            df_biosample = df_biosample.merge(df3, left_on='variant_id', right_index=True)
            # TODO Verify this calculation. Where does this 3 comes from?
            do_not_pass_variant_id_list = df_biosample.ix[df_biosample.variant_count < ((1 / 3) * replicate_count)].variant_id.unique().tolist()
            #
            self.passed_variant_df.loc[do_not_pass_variant_id_list, 'passed'] = False
            self.passed_variant_df.loc[do_not_pass_variant_id_list, this_filter_name] = False
        logger.debug(
            "file: {}; line: {}; Nb variants passed {}".format(__file__, inspect.currentframe().f_lineno,
                                                           (self.passed_variant_df.sum(axis=1) == self.passed_variant_df.shape[1]).sum()))


    def f9_pcr_error(self, var_prop, pcr_error_by_sample):
        """
        Function used to eliminate variants from the data frame in which a pcr error is spotted

        :param var_prop: float, in the [0, 1] interval. It is the threshold chosen by the user
        :param pcr_error_by_sample: boolean, default True. It is used to choose if the analysis is made by sample or sample_replicate
        :return: None
        """
        if not pcr_error_by_sample:
            # TODO: Must be updated for new FilterRunner class
            return
        this_filter_name = inspect.stack()[0][3]
        logger.debug(
            "file: {}; line: {}; {}".format(__file__, inspect.currentframe().f_lineno, this_filter_name))
        if pcr_error_by_sample:
            sample_list = self.variant_read_count_df.biosample_id.unique().tolist()
        else:
            sample_list = None

        for element in sample_list:
            if pcr_error_by_sample:
                df_biosample = self.variant_read_count_df.loc[self.variant_read_count_df['biosample_id'] == element]
                variant_df_by_biosample = df_biosample.merge(self.variant_df, left_on='variant_id', right_on='id')[
                    ['id', 'sequence']].drop_duplicates()
            else:
                df_biosample = self.variant_read_count_df.loc[self.variant_read_count_df['sample_replicate'] == element]
                subset_fasta = os.path.join(self.tempdir, ('{}_{}.fasta'.format(element)))
                variant_id_to_sequence_df = None
            if variant_df_by_biosample.shape[0] > 0:
                ###################################################################
                # 1. Make a fasta file with all variants of the sample or replicate
                ###################################################################
                PathFinder.mkdir_p(os.path.join(self.tempdir, this_filter_name))
                subset_fasta = os.path.join(self.tempdir, this_filter_name, '{}.fasta'.format(element))
                with open(subset_fasta, 'w') as fout:
                    for row in variant_df_by_biosample.itertuples():
                        variant_id = row[0]
                        variant_sequence = row[1]
                        fout.write(">{}\n{}\n".format(variant_id, variant_sequence))
                ###################################################################
                # 2. Determine id threshold for vsearch
                ###################################################################
                L = int(min([len(variant_sequence) for variant_sequence in variant_df_by_biosample.sequence.tolist()]))
                id = floor(((L-1)/L)*100)/100
                ###################################################################
                # 3 Detect all pairs of variants with only 1 difference in the sequences and strong difference in abundance (readcounts)
                # 3.1 vsearch
                ###################################################################
                # import pdb; pdb.set_trace()
                sample_tsv = os.path.join(self.tempdir, '{}.tsv'.format(element))
                vsearch_usearch_global_args = {'db': subset_fasta,
                                               'usearch_global': subset_fasta,
                                               'id': str(id),
                                               'maxrejects': 0,
                                               'maxaccepts': 0,
                                               'userout': sample_tsv,
                                               'userfields': "query+target+alnlen+ids+mism+gaps",
                                               }
                vsearch_usearch_global = VSearch1(**vsearch_usearch_global_args)
                vsearch_usearch_global.run()
                ###################################################################
                # 3.2 If (mism+gaps) ==1
                ###################################################################
                column_names = ['query', 'target', 'alnlen', 'ids', 'mism', 'gaps']
                vsearch_output = pandas.read_csv(sample_tsv, sep='\t', names=column_names)
                # import pdb; pdb.set_trace()
                false_df2 = list(vsearch_output.ix[(vsearch_output.mism + vsearch_output.gaps) != 1].index)
                vsearch_output.drop(false_df2, inplace=True)
                df3 = vsearch_output[['query', 'target']]
                df4 = df_biosample[['variant_id', 'read_count']]
                df3 = df3.merge(df4, left_on=['query'], right_on=['variant_id'])
                df3 = df3.merge(df4, left_on=['target'], right_on=['variant_id'])
                df3['read_count_ratio'] = df3.read_count_x / df3.read_count_y
                df3 = df3[['query', 'read_count_ratio']]
                df3.index = df3['query']
                df_biosample = df_biosample.merge(df3, left_on='variant_id', right_index=True)
                do_not_pass_variant_id_list = list(df_biosample.ix[df_biosample.read_count_ratio < var_prop].variant_id.unique())
                #
                self.passed_variant_df.loc[do_not_pass_variant_id_list, 'passed'] = False
                self.passed_variant_df.loc[do_not_pass_variant_id_list, this_filter_name] = False
        logger.debug(
            "file: {}; line: {}; Nb variants passed {}".format(__file__, inspect.currentframe().f_lineno,
                                                               (self.passed_variant_df.sum(axis=1) ==
                                                                self.passed_variant_df.shape[1]).sum()))

    def f10_chimera(self, chimera_by_sample_replicate, engine, variant_model):
        """
        Function using Vsearch to identify f10_chimera or borderline variant

        :param replicate_obj_list: list of the sample and the corresponding sample replicate
        :param variant_read_count_df: data frame containing the information
        :param chimera_by_sample_replicate: parameter to choose between sample, sample replicate or all
        :return: None
        """
        # TODO: Variant sequence must be taken from variant_df instead of from engine and variant_model
        this_filter_name = inspect.stack()[0][3]
        logger.debug(
            "file: {}; line: {}; {}".format(__file__, inspect.currentframe().f_lineno, this_filter_name))
        variant_id_list_of_lists = []
        if chimera_by_sample_replicate:
            for row in self.variant_read_count_df[['biosample_id', 'replicate_id']].drop_duplicates(['biosample_id', 'replicate_id']).itertuples():
                biosample_id = row.biosample_id
                replicate_id = row.replicate_id
                subset_name = "{}_{}".format(biosample_id, replicate_id)
                df_subset = self.variant_read_count_df.loc[(self.variant_read_count_df.biosample_id == biosample_id) & (self.variant_read_count_df.replicate_id == replicate_id)]
                variant_id_list_of_lists.append(df_subset.variant_id.unique().tolist())
        else:
            for row in self.variant_read_count_df[['biosample_id']].drop_duplicates(['biosample_id']).itertuples():
                biosample_id = row.biosample_id
                subset_name = "{}".format(biosample_id)
                df_subset = self.variant_read_count_df.loc[(self.variant_read_count_df.biosample_id == biosample_id)]
                variant_id_list_of_lists.append(df_subset.variant_id.unique().tolist())
        for variant_id_list in variant_id_list_of_lists:
            if len(variant_id_list) > 0:
                ###################################################################
                # 1. Make a fasta file with all variants of the sample or replicate
                ###################################################################
                PathFinder.mkdir_p(os.path.join(self.tempdir, this_filter_name))
                subset_fasta = os.path.join(self.tempdir, this_filter_name, '{}.fasta'.format(subset_name))
                subset_sortbysize_fasta = os.path.join(self.tempdir, this_filter_name, '{}_sortbysize.fasta'.format(subset_name))
                ###################################################################
                # 2. Sort variants by abundance
                ###################################################################
                with open(subset_fasta, 'w') as fout:
                    variant_sequence_list = []
                    for variant_id in variant_id_list:
                        with engine.connect() as conn:
                            stmt = select([variant_model.__table__.c.sequence]).where(
                                variant_model.__table__.c.id == variant_id)
                            variant_sequence = conn.execute(stmt).first()[0]
                            variant_sequence_list.append(variant_sequence)
                            fout.write(">{}\n{}\n".format(variant_id, variant_sequence))
                vsearch_sortbysize_args = {"sortbysize": subset_fasta, "output": subset_sortbysize_fasta}
                vsearch_sortbysize = Vsearch2(**vsearch_sortbysize_args)
                vsearch_sortbysize.run()
                ###################################################################
                # 3. Run uchime_denovo
                ###################################################################
                subset_borderline_fasta = os.path.join(self.tempdir, this_filter_name,
                                                       '{}_borderline.fasta'.format(subset_name))
                subset_nonchimeras_fasta = os.path.join(self.tempdir, this_filter_name,
                                                       '{}_nonchimeras.fasta'.format(subset_name))
                subset_chimeras_fasta = os.path.join(self.tempdir, this_filter_name,
                                                       '{}_chimeras.fasta'.format(subset_name))
                vsearch_chimera_args = {
                    "uchime_denovo": subset_sortbysize_fasta,
                    "borderline": subset_borderline_fasta,
                    "nonchimeras": subset_nonchimeras_fasta,
                    "chimeras": subset_chimeras_fasta
                }
                vsearch_chimera = Vsearch3(**vsearch_chimera_args)
                vsearch_chimera.run()
                ###################################################################
                # 4. Delete variant from replicate/sample if chimeras
                ###################################################################
                subset_chimera_seqio = SeqIO.parse(open(subset_chimeras_fasta), 'fasta')
                for subset_chimera in subset_chimera_seqio:
                    self.passed_variant_df.loc[variant_id, 'passed'] = False
                    self.passed_variant_df.loc[subset_chimera.id, 'f10_chimera'] = False # does not pass
                ###################################################################
                # 5. Flag borderline sequences as possible chimera
                ###################################################################
                subset_borderline_seqio = SeqIO.parse(open(subset_borderline_fasta), 'fasta')
                for subset_borderline in subset_borderline_seqio:
                    self.passed_variant_df.loc[subset_borderline.id, 'f10_pass_chimera_borderline'] = False # does not pass
        logger.debug(
            "file: {}; line: {}; Nb variants passed {}".format(__file__, inspect.currentframe().f_lineno,
                                                               (self.passed_variant_df.sum(axis=1) ==
                                                                self.passed_variant_df.shape[1]).sum()))


    def f11_renkonen(self, number_of_replicate, renkonen_tail):
        """
        Function calculating distance between sample replicate of the sample to determine if the sample replicate
        must be deleted or not
        :param engine: engine of the database
        :param replicate_model: model of the replicate table
        :param variant_read_count_df: data frame containing the information
        :param number_of_replicate: Number of replicate by sample
        :return: None
        """
        # TODO: Must be updated for new FilterRunner class
        return
        this_filter_name = inspect.stack()[0][3]
        logger.debug(
            "file: {}; line: {}; {}".format(__file__, inspect.currentframe().f_lineno, this_filter_name))
        ########################################################
        # proportion of the reads of variant i per replicate j (Ni,j=1/Nj=1)
        ########################################################
        variant_read_proportion_per_replicate_df = self.variant_read_count_df[['biosample_id', 'replicate_id', 'read_count']].groupby(
            by=['biosample_id', 'replicate_id']).sum().reset_index()
        # Merge the column with the total reads by sample replicates for calculate the ratio
        variant_read_proportion_per_replicate_df = self.variant_read_count_df.merge(variant_read_proportion_per_replicate_df, left_on=['biosample_id', 'replicate_id'],
                                               right_on=['biosample_id', 'replicate_id'])
        variant_read_proportion_per_replicate_df.columns = ['variant_id', 'biosample_id', 'replicate_id',
                       'read_count_per_variant_per_biosample_per_replicate', 'read_count_per_biosample_replicate']
        variant_read_proportion_per_replicate_df[
            'read_proportion_of_variant_in_replicate'] = variant_read_proportion_per_replicate_df.read_count_per_variant_per_biosample_per_replicate / variant_read_proportion_per_replicate_df.read_count_per_biosample_replicate
        # for biosample_id in self.variant_read_count_df.biosample_id.unique():
        # replicate_combinatorics = itertools.permutations(self.variant_read_count_df.replicate_id.unique().tolist(), 2)
        biosample_id = 1
        replicate_id1 = 1
        replicate_id2 = 2
        variant_read_proportion_per_replicate_per_biosample_df = variant_read_proportion_per_replicate_df.loc[
            variant_read_proportion_per_replicate_df.biosample_id == biosample_id]
        ########################################################
        # 2. Calculate renkonen distance index (D) for all pairs of replicates of the same sample
        ########################################################
        variant_read_proportion_per_replicate1_per_biosample_df = variant_read_proportion_per_replicate_per_biosample_df.loc[
            variant_read_proportion_per_replicate_per_biosample_df.replicate_id == replicate_id1, ['variant_id', 'replicate_id', 'read_proportion_of_variant_in_replicate']]
        variant_read_proportion_per_replicate2_per_biosample_df = variant_read_proportion_per_replicate_per_biosample_df.loc[
            variant_read_proportion_per_replicate_per_biosample_df.replicate_id == replicate_id2, ['variant_id', 'replicate_id', 'read_proportion_of_variant_in_replicate']]
        variant_read_proportion_per_replicate_1_2 = variant_read_proportion_per_replicate1_per_biosample_df.merge(
            variant_read_proportion_per_replicate2_per_biosample_df, on='variant_id')
        variant_read_proportion_per_replicate_1_2.columns = ['variant_id', 'replicate_id1',
                                                             'read_proportion_of_variant_in_replicate1',
                                                             'replicate_id2',
                                                             'read_proportion_of_variant_in_replicate_2']
        variant_read_proportion_per_replicate_1_2 = variant_read_proportion_per_replicate_1_2[
            ['variant_id', 'read_proportion_of_variant_in_replicate1', 'read_proportion_of_variant_in_replicate2']]
        variant_read_proportion_per_replicate_1_2['min_read_proportion'] = variant_read_proportion_per_replicate_1_2[
            ['read_proportion_of_variant_in_replicate1', 'read_proportion_of_variant_in_replicate2']].apply(
            lambda row: row.min(), axis=1)
        #
        columns_name = ['repl_i', 'repl_j', 'distance']
        df_read_count_per_sample_replicate = self.variant_read_count_df.groupby(by=['sample_replicate'])['count'].sum()
        df_read_count_per_sample_replicate = df_read_count_per_sample_replicate.to_frame()
        df_read_count_per_sample_replicate.columns = ['replicate_count']
        df_read_count_per_sample_replicate = self.variant_read_count_df.merge(df_read_count_per_sample_replicate, left_on='sample_replicate', right_index=True)
        df_read_count_per_sample_replicate['proportion'] = df_read_count_per_sample_replicate['count'] / df_read_count_per_sample_replicate['replicate_count']
        # df_replicate = df_read_count_per_sample_replicate.groupby(by=['biosample'])['sample_replicate'].to_frame()
        samples = df_read_count_per_sample_replicate['biosample']
        samples = list(set(samples.tolist()))
        for sample in samples:
            df_permutation_distance = pandas.DataFrame(columns=columns_name)
            df_replicate = df_read_count_per_sample_replicate.loc[df_read_count_per_sample_replicate['biosample'] == sample]
            replicates = list(set(df_replicate['sample_replicate'].tolist()))
            for combi in itertools.permutations(replicates, 2):
                combi = list(combi)
                df_repli = df_replicate.loc[df_replicate['sample_replicate'] == combi[0]]
                data_repli = df_repli[['variant_seq', 'sample_replicate', 'proportion']]
                df_replj = df_replicate.loc[df_replicate['sample_replicate'] == combi[1]]
                data_replj = df_replj[['variant_seq', 'sample_replicate', 'proportion']]
                df_replij = data_repli.append(data_replj)
                group_repl = df_replij.groupby(by=['variant_seq'])['proportion'].min()
                distance = 1 - group_repl.sum()
                query = [combi[0], combi[1], distance]
                df_permutation_distance.loc[len(df_permutation_distance)] = query
            # df_calc = df_permutation_distance.loc[df_permutation_distance['repl_i'] == combi[0]]
            indices_to_drop = list(
                df_permutation_distance.loc[df_permutation_distance.distance > renkonen_tail].index.tolist()
            )
            df_permutation_distance.drop(indices_to_drop, inplace=True)
            repl_list = list(set(df_permutation_distance['repl_i'].tolist()))
            for repl_i in repl_list:
                df_calc = df_permutation_distance.loc[df_permutation_distance['repl_i'] == repl_i]
                if len(df_calc) > ((number_of_replicate -1)/2):
                    index_to_drop = self.variant_read_count_df.loc[self.variant_read_count_df['sample_replicate'] == repl_i].index.tolist()
                    self.passsed_variant_ids = sorted(list(set(index_to_drop + self.passsed_variant_ids)))


    def f12_indel(self):
        """

        :param delete_var:
        :return:
        """
        this_filter_name = inspect.stack()[0][3]
        logger.debug(
            "file: {}; line: {}; {}".format(__file__, inspect.currentframe().f_lineno, this_filter_name))
        df = self.variant_df.copy()
        df['sequence_length_module_3'] = self.variant_df.sequence.apply(lambda x: len(x) % 3)
        majority_sequence_length_module_3 = df.sequence_length_module_3.mode()
        df = df.loc[df['sequence_length_module_3'] != majority_sequence_length_module_3.values[0]]
        do_not_pass_variant_id_list = df.id.tolist()
        #
        # self.passed_variant_df.loc[do_not_pass_variant_id_list, 'passed'] = False
        self.passed_variant_df.loc[do_not_pass_variant_id_list, this_filter_name] = False
        logger.debug(
            "file: {}; line: {}; Nb variants passed {}".format(__file__, inspect.currentframe().f_lineno,
                                                               (self.passed_variant_df.sum(axis=1) == self.passed_variant_df.shape[1]).sum()))


    def codon_stop(self, df_codon_stop_per_genetic_code, genetic_code, delete_var):
        """
        Function searching stop codon the different reading frame of a sequence and tag/delete variant if there is
        stop codon in the 3 reading frames
        :param df_codon_stop_per_genetic_code: data frame which contains all codon stop classified by genetic code
        :param genetic_code: genetic code to search stop codon in the df_codon_stop_per_genetic_code dataframe
        :param delete_var: option which define if the variants must be deleted or tagged
        :return: void
        """
        # TODO: Must be reviewed for new object
        return
        sequences = self.variant_read_count_df['variant_seq'].tolist()
        # Get the list of codon stop according to the genetic code given by the user
        df2 = df_codon_stop_per_genetic_code.loc[df_codon_stop_per_genetic_code['genetic_code'] == genetic_code]
        codon_stop_list = df2['codon'].tolist()
        indices_to_drop = []
        # Take sequence 1 by 1
        for sequence in sequences:
            # Keep the entire sequence for select in the dataframe
            entire_sequence = sequence
            # define the current reading_frame
            reading_frame = 1
            # Number of times when a codon stop is fin in the sequence
            fail = 0
            for i in range(3):
                # if a codon stop is not find in the first reading frame no need to search the others
                if reading_frame > 1 and fail == 0:
                    break
                # Get a new sequence eliminating the first nucleotide
                if reading_frame == 2:
                    sequence = sequence[1:]
                # Get a new sequence eliminating the first and the second nucleotides
                elif reading_frame == 3:
                    sequence = sequence[2:]
                # Divide the sequence in length 3 codons
                codon_list = re.findall('...', sequence)
                codon_stop_count = 0
                # Count the number of codon stop in the sequence length
                for codon_stop in codon_stop_list:
                    if codon_stop in codon_list:
                        codon_stop_count += 1
                # If a codon stop or more are found then 1 fail is added
                if codon_stop_count > 0:
                    fail += 1
                reading_frame += 1
            # If a codon stop is find in every reading frames, all the indexes where the variant sequence appeared are
            # kept and will be tagged or deleted at the loop end
            if fail == 3:
                indices_to_drop.append(self.variant_read_count_df.loc[self.variant_read_count_df['variant_seq'] == entire_sequence].index().tolist())
        # Depending on the user choice the variant will be tagged in dataframe or removed from it
        if delete_var:
            self.passsed_variant_ids = sorted(list(set(indices_to_drop + self.passsed_variant_ids)))
        else:
            self.variant_read_count_df['is_pseudogene_codon_stop'] = (self.variant_read_count_df.index in indices_to_drop)

    def consensus(self):
        """
        Function used to display the read average of the remaining variant
        :return:
        """
        # TODO: Must be reviewed for new object
        variants_sequences = self.variant_read_count_df['variant_seq']
        variants_sequences =list(set(variants_sequences))
        read_average_columns = ['variant', 'read_average']
        read_average_df = pandas.DataFrame(columns=read_average_columns)
        for variant in variants_sequences:
            variant_df = self.variant_read_count_df.loc[self.variant_read_count_df['variant_seq'] == variant]
            read_average = round(variant_df["count"].sum()/len(variant_df['count']), 0)
            read_average_df.loc[len(read_average_df)] = [variant, read_average]
        self.variant_read_count_df = self.variant_read_count_df.merge(read_average_df, left_on='variant_seq', right_on='variant')
        self.variant_read_count_df = self.variant_read_count_df.drop(columns=['variant'])

    def filtered_variants(self):
        return self.variant_read_count_df
