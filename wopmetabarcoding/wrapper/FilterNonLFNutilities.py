# -*- coding: utf-8 -*-
"""NON LFN Filters

This module will store and run these NONLFNFilters:

"""

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





class FilterNonLFNRunner:


    def __init__(self, variant_read_count_df,marker_id,variant_df):
        logger.debug(
            "file: {}; line: {}; FilterNonLFNRunner.__init__".format(__file__, inspect.currentframe().f_lineno))
        # self.variant_df = variant_df
        self.variant_read_count_df = variant_read_count_df[['marker_id', 'run_id', 'variant_id', 'biosample_id', 'replicate_id', 'read_count']]
        self.variant_df = variant_df[['id','sequence']]
        self.marker_id = marker_id
        #
        self.tempdir = os.path.join(tempdir, "FilterUtilities", "FilterUtilities", self.__class__.__name__)
        PathFinder.mkdir_p(self.tempdir)
        #
        if self.variant_read_count_df.shape[1] != 6:
            raise Exception('Columns missing in the variant2sample2replicate2count data frame!')

        # if self.variant_df.shape[1] != 2:
        #     raise Exception('Columns missing in the variant2sample2replicate2count data frame!')
        #
        ################################
        # Output df with deleted variants
        ################################
        self.delete_variant_df = pandas.DataFrame(data={'run_id':[],'marker_id':[],'variant_id':[], 'biosample_id':[], 'replicate_id':[], 'read_count':[]}, dtype='int')
        logger.debug(
            "file: {}; line: {}; Initial nb of variants {}".format(__file__, inspect.currentframe().f_lineno,
                                                               (self.delete_variant_df.sum(axis=1) == self.delete_variant_df.shape[1]).sum()))




    def f10_pcr_error(self, pcr_error_var_prop):
        """
        Function used to eliminate variants from the data frame in which a pcr error is spotted

        :param var_prop: float, in the [0, 1] interval. It is the threshold chosen by the user
        :param pcr_error_by_sample: boolean, default True. It is used to choose if the analysis is made by sample or sample_replicate
        :return: None
        """
        this_filter_name ='f10_pcr_error'
        ###################################################################
        # 1. Make a fasta file with all variants of the sample or replicate
        ###################################################################
        #
        # PathFinder.mkdir_p(os.path.join(self.tempdir, this_filter_name))
        # variant_fasta = os.path.join(self.tempdir, this_filter_name, '{}.fasta'.format("pcr_error"))
        # with open(variant_fasta, 'w') as fout:
        #     for row in self.variant_df.itertuples():
        #         id = row[0]
        #         sequence = row[1]
        #         fout.write(">{}\n{}\n".format(id, sequence))

        PathFinder.mkdir_p(os.path.join(self.tempdir, this_filter_name))
        variant_fasta = os.path.join(self.tempdir, this_filter_name, '{}.fasta'.format("pcr_error"))
        with open(variant_fasta, 'w') as fout:
            for row in self.variant_df.itertuples():
                variant_id = row[0]
                variant_sequence = row[1]
                fout.write(">{}\n{}\n".format(variant_id, variant_sequence))

        ###################################################################
        # 2. Determine id threshold for vsearch
        ###################################################################
        L = int(min([len(sequence) for sequence in self.variant_df.sequence.tolist()]))
        id = floor(((L-1)/L)*100)/100
        ###################################################################
        # 3 Detect all pairs of variants with only 1 difference in the sequences and strong difference in abundance (readcounts)
        # 3.1 vsearch
        ###################################################################
        # import pdb; pdb.set_trace()
        sample_tsv = os.path.join(self.tempdir, '{}.tsv'.format("pcr_error"))
        vsearch_usearch_global_args = {'db': variant_fasta,
                                       'usearch_global': variant_fasta,
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
        df_read_count = self.variant_read_count_df
        df4 = df_read_count[['variant_id', 'read_count']]
        df3 = df3.merge(df4, left_on=['query'], right_on=['variant_id'])
        df3 = df3.merge(df4, left_on=['target'], right_on=['variant_id'])
        df3['read_count_ratio'] = df3.read_count_x / df3.read_count_y
        df3 = df3[['query', 'read_count_ratio']]
        df3.index = df3['query']
        list_variant_delete = list(df3.loc[df3.read_count_ratio < pcr_error_var_prop].variant_id_x.unique())
        passed_pcr_error_df = self.variant_read_count_df.loc[self.variant_read_count_df.variant_id != list]

        #  Prepare output df and concatenate to self.delete_variant_df
        self.delete_variant_df = pandas.concat([self.delete_variant_df, passed_pcr_error_df], sort=False)