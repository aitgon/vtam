from itertools import tee
from unittest import TestCase

import itertools

from wopmetabarcoding.utils.PathFinder import PathFinder
import os
from wopmetabarcoding.utils.constants import tempdir
import pandas
from wopmetabarcoding.wrapper.FilterRenkonen import renkonen_distance, renkonen_distance_df


class TestFilterRenkonen(TestCase):

    def setUp(self):
        self.variant_read_count_df = pandas.DataFrame({
            'run_id': [1] * 6,
            'marker_id': [1] * 6,
            'variant_id': [1] * 3 + [2] * 3,
            'biosample_id': [1] * 6,
            'replicate_id': [1, 2, 3] * 2,
            'read_count': [
                5180, 5254, 9378, 193, 99, 209
            ],
        })

        self.tempdir = os.path.join(tempdir, "FilterUtilities", self.__class__.__name__)
        PathFinder.mkdir_p(self.tempdir)

    def test_f12_delete_filter_renkonen(self):
        #
        Rthr = 0.005
        #
        df2 = self.variant_read_count_df.groupby(['run_id', 'marker_id', 'biosample_id']).agg('replicate_id').apply(lambda x: list(set(x))).reset_index()
        df2['threshold_distance_number'] = df2['replicate_id'].apply(lambda x: (len(x) - 1) / 2 )
        df2['replicate_id_pairwise'] = df2.replicate_id.apply(lambda x: list(itertools.combinations(x, 2)))
        df2.drop('replicate_id', axis=1, inplace=True)
        df3 = pandas.DataFrame(data={'run_id': [], 'marker_id': [], 'biosample_id': [], 'left_replicate_id': [], 'right_replicate_id': [], 'renkonen_distance': []}, dtype='int')
        for row in df2.itertuples():
            run_id = row.run_id
            marker_id = row.marker_id
            biosample_id = row.biosample_id
            replicate_id_pairwise = row.replicate_id_pairwise
            for left_replicate_id, right_replicate_id in replicate_id_pairwise:
                df3 = pandas.concat([df3, pandas.DataFrame({'run_id': [run_id], 'marker_id': [marker_id], 'biosample_id': [biosample_id],
                            'left_replicate_id': [left_replicate_id], 'right_replicate_id': [right_replicate_id]})], axis=0)
        #
        for row in df3.itertuples():
            run_id = row.run_id
            marker_id = row.marker_id
            biosample_id = row.biosample_id
            left_replicate_id = row.left_replicate_id
            right_replicate_id = row.right_replicate_id
            D = renkonen_distance(self.variant_read_count_df, run_id, marker_id, biosample_id, left_replicate_id,
                              right_replicate_id)
            df3.loc[(df3.run_id==run_id) & (df3.marker_id==marker_id) & (df3.biosample_id==biosample_id)
                & (df3.left_replicate_id==left_replicate_id) & (df3.right_replicate_id==right_replicate_id),'renkonen_distance'] = D
        #
        df3['is_distance_gt_rthr']=df3.renkonen_distance > Rthr
        #
        df4 = pandas.DataFrame(
            data={'run_id': [], 'marker_id': [], 'biosample_id': [], 'replicate_id': [], 'is_distance_gt_rthr': []}, dtype='int')
        df4 = df4.rename(columns={'replicate_id': 'left_replicate_id'})
        df4 = pandas.concat([df4, df3[['run_id', 'marker_id', 'biosample_id', 'left_replicate_id','is_distance_gt_rthr']]])
        df4 = df4.rename(columns={'left_replicate_id': 'right_replicate_id'})
        df4 = pandas.concat([df4, df3[['run_id', 'marker_id', 'biosample_id', 'right_replicate_id','is_distance_gt_rthr']]], axis=0)
        df4 = df4.rename(columns={'right_replicate_id': 'replicate_id'})
        #
        df5 = df4.groupby(['run_id', 'marker_id', 'biosample_id', 'replicate_id']).sum().reset_index()
        df5 = df5.rename(columns={'is_distance_gt_rthr': 'distance_number'})
        df5 = df5.merge(df2[['run_id', 'marker_id', 'biosample_id', 'threshold_distance_number']])
        #
        df5['filter_delete'] = False
        df5.loc[df5.distance_number > df5.threshold_distance_number, 'filter_delete'] = True
        #
        dfout = self.variant_read_count_df.merge(df5)
        dfout.drop(['distance_number', 'threshold_distance_number'], axis=1, inplace=True)
        dfout['filter_id'] = 12


    def test_renkonen_distance(self):
        # Output
        # biosample_id 1, replicate_id 1, replicate_id 2, renkonen_similarity and distance 0.982573959807409 and 0.017426040192591
        # biosample_id 1, replicate_id 1, replicate_id 3, renkonen_similarity and distance 0.985880012193912 and 0.014119987806088
        run_id = 1
        marker_id = 1
        biosample_id = 1
        left_replicate_id = 1
        right_replicate_id = 2
        #
        distance_left_right = renkonen_distance(self.variant_read_count_df,run_id,marker_id,biosample_id,left_replicate_id,right_replicate_id)
        self.assertAlmostEqual(distance_left_right, 0.017426040192591)
        right_replicate_id = 3
        distance_left_right= renkonen_distance(self.variant_read_count_df,run_id,marker_id,biosample_id,left_replicate_id,right_replicate_id)
        self.assertAlmostEqual(distance_left_right, 0.014119987806088)



    def test_delete_replicate(self):
        # TODO

        # Input
        variant_read_count_df = pandas.DataFrame({
            'run_id': [1] * 6,
            'marker_id': [1] * 6,
            'variant_id': [1] * 3 + [2] * 3,
            'biosample_id': [1] * 6,
            'replicate_id': [1, 2, 3] * 2,
            'read_count': [
                5180, 5254, 9378, 193, 99, 209
            ],
        })






    def test_f11_filter_renkonen(self):
            """

            :return: None
            """
            #

            # logger.debug(
            #     "file: {}; line: {}; {}".format(__file__, inspect.currentframe().f_lineno, this_filter_name))


            ########################################################
            # proportion of the reads of variant i per replicate j (Ni,j=1/Nj=1)
            ########################################################

            renkonen_tail = 0.013
            number_of_replicate = 2
            passsed_variant_ids= []
            variant_read_proportion_per_replicate_df = self.variant_read_count_df[['run_id', 'marker_id', 'biosample_id', 'replicate_id', 'read_count']].groupby(
                                                                                            by=['run_id', 'marker_id', 'biosample_id', 'replicate_id']).sum().reset_index()

            # # Merge the column with the total reads by sample replicates for calculate the ratio
            # variant_read_proportion_per_replicate_df = self.variant_read_count_df.merge(variant_read_proportion_per_replicate_df, left_on=['biosample_id', 'replicate_id'],
            #                                                                             right_on=['biosample_id', 'replicate_id'])
            # variant_read_proportion_per_replicate_df.columns = ['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id','rc_per_v_per_b_per_r',
            #                                                        'rc_per_b_r']
            #
            #
            # variant_read_proportion_per_replicate_df['rp_of_variant_in_replicate'] = variant_read_proportion_per_replicate_df.rc_per_v_per_b_per_r / variant_read_proportion_per_replicate_df.rc_per_b_r
            #
            # # for biosample_id in self.variant_read_count_df.biosample_id.unique():
            # # replicate_combinatorics = itertools.permutations(self.variant_read_count_df.replicate_id.unique().tolist(), 2)
            # biosample_id = 1
            # replicate_id1 = 1
            # replicate_id2 = 2
            # variant_read_proportion_per_replicate_per_biosample_df = variant_read_proportion_per_replicate_df.loc[
            #                                                                 variant_read_proportion_per_replicate_df.biosample_id == biosample_id]
            #
            # ########################################################
            # # 2. Calculate renkonen distance index (D) for all pairs of replicates of the same sample
            # ########################################################
            # variant_read_proportion_per_replicate1_per_biosample_df = variant_read_proportion_per_replicate_per_biosample_df.loc[
            #     variant_read_proportion_per_replicate_per_biosample_df.replicate_id == replicate_id1, ['variant_id', 'replicate_id', 'rp_of_variant_in_replicate']]
            # variant_read_proportion_per_replicate2_per_biosample_df = variant_read_proportion_per_replicate_per_biosample_df.loc[
            #     variant_read_proportion_per_replicate_per_biosample_df.replicate_id == replicate_id2, ['variant_id', 'replicate_id', 'rp_of_variant_in_replicate']]
            # variant_read_proportion_per_replicate_1_2 = variant_read_proportion_per_replicate1_per_biosample_df.merge(variant_read_proportion_per_replicate2_per_biosample_df,
            #                                                                                                           on='variant_id')
            # variant_read_proportion_per_replicate_1_2.columns = ['variant_id', 'replicate_id1',
            #                                                      'rp_of_variant_in_replicate1',
            #                                                      'replicate_id2',
            #                                                      'rp_of_variant_in_replicate_2']
            #
            # variant_read_proportion_per_replicate_1_2 = variant_read_proportion_per_replicate_1_2[['variant_id','rp_of_variant_in_replicate1', 'rp_of_variant_in_replicate_2']]
            #
            #
            # variant_read_proportion_per_replicate_1_2['min_read_proportion'] = variant_read_proportion_per_replicate_1_2[
            #     ['rp_of_variant_in_replicate1', 'rp_of_variant_in_replicate_2']].apply(lambda row: row.min(), axis=1)
            #
            # columns_name = ['repl_i', 'repl_j', 'distance']
            # df_read_count_per_sample_replicate = self.variant_read_count_df.groupby(by=['replicate_id'])['read_count'].sum()
            # df_read_count_per_sample_replicate = df_read_count_per_sample_replicate.to_frame()
            # df_read_count_per_sample_replicate.columns = ['replicate_count']
            # df_read_count_per_sample_replicate = self.variant_read_count_df.merge(df_read_count_per_sample_replicate, left_on='replicate_id', right_index=True)
            # df_read_count_per_sample_replicate['proportion'] = df_read_count_per_sample_replicate['read_count'] / df_read_count_per_sample_replicate['replicate_count']
            #
            # # df_replicate = df_read_count_per_sample_replicate.groupby(by=['biosample'])['sample_replicate'].to_frame()
            # samples = df_read_count_per_sample_replicate['biosample_id']
            # samples = list(set(samples.tolist()))
            # #
            # for sample in samples:
            #     df_permutation_distance = pandas.DataFrame(columns=columns_name)
            #     df_replicate = df_read_count_per_sample_replicate.loc[df_read_count_per_sample_replicate['biosample_id'] == sample]
            #     replicates = list(set(df_replicate['replicate_id'].tolist()))
            #     for combi in itertools.permutations(replicates, 2):
            #         combi = list(combi)
            #         df_repli = df_replicate.loc[df_replicate['replicate_id'] == combi[0]]
            #         # import pdb;
            #         # pdb.set_trace()
            #         data_repli = df_repli[['variant_id', 'replicate_id', 'proportion']]
            #         df_replj = df_replicate.loc[df_replicate['replicate_id'] == combi[1]]
            #         data_replj = df_replj[['variant_id', 'replicate_id', 'proportion']]
            #         df_replij = data_repli.append(data_replj)
            #         group_repl = df_replij.groupby(by=['variant_id'])['proportion'].min()
            #         distance = 1 - group_repl.sum()
            #         query = [combi[0], combi[1], distance]
            #         df_permutation_distance.loc[len(df_permutation_distance)] = query
            #     # df_calc = df_permutation_distance.loc[df_permutation_distance['repl_i'] == combi[0]]
            #     indices_to_drop = list(
            #         df_permutation_distance.loc[df_permutation_distance.distance > renkonen_tail].index.tolist()
            #     )
            #     df_permutation_distance.drop(indices_to_drop, inplace=True)
            #     repl_list = list(set(df_permutation_distance['repl_i'].tolist()))
            #     for repl_i in repl_list:
            #         df_calc = df_permutation_distance.loc[df_permutation_distance['repl_i'] == repl_i]
            #         if len(df_calc) > ((number_of_replicate -1) / 2):
            #             index = self.variant_read_count_df.loc[self.variant_read_count_df['replicate_id'] == repl_i].index.tolist()
            #             passsed_variant_ids = sorted(list(set(index + passsed_variant_ids)))
            # # import pdb;
            # # pdb.set_trace()
            # return passsed_variant_ids