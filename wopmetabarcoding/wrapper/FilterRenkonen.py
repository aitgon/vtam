import time

import pandas, itertools
from sqlalchemy import select
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper

from wopmetabarcoding.utils.utilities import create_step_tmp_dir


class FilterRenkonen(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.FilterRenkonen"
    }

    # Input file
    __input_file_sample2fasta = "sample2fasta"
    # Input table
    __input_table_marker = "Marker"
    __input_table_run = "Run"
    __input_table_biosample = "Biosample"
    __input_table_replicate = "Replicate"
    __input_table_chimera = "FilterChimera"
    # Output table
    __output_table_filter_renkonen = "FilterRenkonen"

    def specify_input_file(self):
        return [
            FilterRenkonen.__input_file_sample2fasta,

        ]

    def specify_input_table(self):
        return [
            FilterRenkonen.__input_table_marker,
            FilterRenkonen.__input_table_run,
            FilterRenkonen.__input_table_biosample,
            FilterRenkonen.__input_table_replicate,
            FilterRenkonen.__input_table_chimera,
        ]

    def specify_output_table(self):
        return [
            FilterRenkonen.__output_table_filter_renkonen,
        ]

    def specify_params(self):
        return {
            "renkonen_threshold": "float",
        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind

        ##########################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ##########################################################
        #
        # Input file path
        input_file_sample2fasta = self.input_file(FilterRenkonen.__input_file_sample2fasta)
        #
        # Input table models
        marker_model = self.input_table(FilterRenkonen.__input_table_marker)
        run_model = self.input_table(FilterRenkonen.__input_table_run)
        chimera_model = self.input_table(FilterRenkonen.__input_table_chimera)
        biosample_model = self.input_table(FilterRenkonen.__input_table_biosample)
        replicate_model = self.input_table(FilterRenkonen.__input_table_replicate)
        #
        # Options
        # TaxAssign_bak parameters
        renkonen_threshold = float(self.option("renkonen_threshold"))
        #
        # Output table models
        filter_renkonen_model = self.output_table(FilterRenkonen.__output_table_filter_renkonen)

        ##########################################################
        #
        # 1. Read sample2fasta to get run_id, marker_id, biosample_id, replicate_id for current analysis
        #
        ##########################################################
        sample2fasta_df = pandas.read_csv(input_file_sample2fasta, sep="\t", header=None, \
                                          names=['tag_forward', 'primer_forward', 'tag_reverse', 'primer_reverse',
                                                 'marker_name', 'biosample_name', \
                                                 'replicate_name', 'run_name', 'fastq_fwd', 'fastq_rev', 'fasta'])
        sample_instance_list = []
        for row in sample2fasta_df.itertuples():
            marker_name = row.marker_name
            run_name = row.run_name
            biosample_name = row.biosample_name
            replicate_name = row.replicate_name
            with engine.connect() as conn:
                # get run_id ###########
                stmt_select_run_id = select([run_model.__table__.c.id]).where(run_model.__table__.c.name == run_name)
                run_id = conn.execute(stmt_select_run_id).first()[0]
                # get marker_id ###########
                stmt_select_marker_id = select([marker_model.__table__.c.id]).where(
                    marker_model.__table__.c.name == marker_name)
                marker_id = conn.execute(stmt_select_marker_id).first()[0]
                # get biosample_id ###########
                stmt_select_biosample_id = select([biosample_model.__table__.c.id]).where(
                    biosample_model.__table__.c.name == biosample_name)
                biosample_id = conn.execute(stmt_select_biosample_id).first()[0]
                # get replicate_id ###########
                stmt_select_replicate_id = select([replicate_model.__table__.c.id]).where(
                    replicate_model.__table__.c.name == replicate_name)
                replicate_id = conn.execute(stmt_select_replicate_id).first()[0]
                # add this sample_instance ###########
                sample_instance_list.append({'run_id': run_id, 'marker_id': marker_id, 'biosample_id': biosample_id,
                                             'replicate_id': replicate_id})

        ##########################################################
        #
        # 2. Delete /run/markerbiosample/replicate from this filter table
        #
        ##########################################################
        with engine.connect() as conn:
            conn.execute(filter_renkonen_model.__table__.delete(), sample_instance_list)

        ##########################################################
        #
        # 3. Select marker/run/biosample/replicate from variant_read_count_model
        #
        ##########################################################

        chimera_model_table = chimera_model.__table__
        stmt_variant_filter_lfn = select([chimera_model_table.c.run_id,
                                          chimera_model_table.c.marker_id,
                                          chimera_model_table.c.biosample_id,
                                          chimera_model_table.c.replicate_id,
                                          chimera_model_table.c.variant_id,
                                          chimera_model_table.c.read_count]) \
            .where(chimera_model_table.c.filter_delete == 0)
        # Select to DataFrame
        variant_filter_lfn_passed_list = []
        with engine.connect() as conn:
            for row in conn.execute(stmt_variant_filter_lfn).fetchall():
                variant_filter_lfn_passed_list.append(row)
        variant_read_count_df = pandas.DataFrame.from_records(variant_filter_lfn_passed_list,
                                                              columns=['run_id', 'marker_id',
                                                                       'biosample_id', 'replicate_id', 'variant_id', 'read_count',])

        ##########################################################
        #
        # 4. Run Filter
        #
        ##########################################################

        df = f12_filter_delete_renkonen(variant_read_count_df, renkonen_threshold)

        ##########################################################
        #
        # 5. Insert Filter data
        #
        ##########################################################
        with engine.connect() as conn:
                conn.execute(filter_renkonen_model.__table__.insert(), df.to_dict('records'))





def  renkonen_distance(variant_read_count_df, run_id, marker_id, biosample_id, left_replicate_id, right_replicate_id):
    #  Compute sum of read_count per 'run_id', 'marker_id', 'biosample_id', 'replicate_id'
    variant_read_proportion_per_replicate_df = variant_read_count_df[
        ['run_id', 'marker_id', 'biosample_id', 'replicate_id', 'read_count']].groupby(
        by=['run_id', 'marker_id', 'biosample_id', 'replicate_id']).sum().reset_index()
    variant_read_proportion_per_replicate_df = variant_read_proportion_per_replicate_df.rename(
        columns={'read_count': 'read_count_sum_per_replicate'})

    # Merge variant read_count with read_count_sum_per_replicate
    variant_read_proportion_per_replicate_df = variant_read_count_df.merge(variant_read_proportion_per_replicate_df,
                                                                           left_on=['run_id', 'marker_id',
                                                                                    'biosample_id', 'replicate_id'],
                                                                           right_on=['run_id', 'marker_id',
                                                                                     'biosample_id', 'replicate_id'])

    variant_read_proportion_per_replicate_df['variant_read_count_propotion_per_replicate'] \
        = variant_read_proportion_per_replicate_df.read_count / variant_read_proportion_per_replicate_df.read_count_sum_per_replicate
    variant_read_proportion_per_replicate_df.drop('read_count', axis=1, inplace=True)
    variant_read_proportion_per_replicate_df.drop('read_count_sum_per_replicate', axis=1, inplace=True)
    #

    # Select the read proportion for the biosample_id, left_replicate
    left_variant_read_proportion_per_replicate_per_biosample_df = variant_read_proportion_per_replicate_df.loc[
        (variant_read_proportion_per_replicate_df.run_id == run_id)
        & (variant_read_proportion_per_replicate_df.marker_id == marker_id)
        & (variant_read_proportion_per_replicate_df.biosample_id == biosample_id)
        & (variant_read_proportion_per_replicate_df.replicate_id == left_replicate_id)
        ]

    # Select the read proportion for the biosample_id, left_replicate
    right_variant_read_proportion_per_replicate_per_biosample_df = variant_read_proportion_per_replicate_df.loc[
        (variant_read_proportion_per_replicate_df.run_id == run_id)
        & (variant_read_proportion_per_replicate_df.marker_id == marker_id)
        & (variant_read_proportion_per_replicate_df.biosample_id == biosample_id)
        & (variant_read_proportion_per_replicate_df.replicate_id == right_replicate_id)
        ]


    #  Merge left and right replicate
    variant_read_proportion_per_replicate_left_right = left_variant_read_proportion_per_replicate_per_biosample_df.merge( \
        right_variant_read_proportion_per_replicate_per_biosample_df,
        on=['run_id', 'marker_id', 'variant_id', 'biosample_id'])
    # rename columns
    variant_read_proportion_per_replicate_left_right = variant_read_proportion_per_replicate_left_right.rename(
        columns={'replicate_id_x': 'replicate_id_left'})
    variant_read_proportion_per_replicate_left_right = variant_read_proportion_per_replicate_left_right.rename(
        columns={'variant_read_count_propotion_per_replicate_x': 'variant_read_count_propotion_per_replicate_left'})
    variant_read_proportion_per_replicate_left_right = variant_read_proportion_per_replicate_left_right.rename(
        columns={'replicate_id_y': 'replicate_id_right'})
    variant_read_proportion_per_replicate_left_right = variant_read_proportion_per_replicate_left_right.rename(
        columns={'variant_read_count_propotion_per_replicate_y': 'variant_read_count_propotion_per_replicate_right'})

    # variant_read_proportion_per_replicate_left_right = variant_read_proportion_per_replicate_left_right[['variant_id','rp_of_variant_in_replicate1', 'rp_of_variant_in_replicate2']]

    variant_read_proportion_per_replicate_left_right['min_read_proportion'] = \
    variant_read_proportion_per_replicate_left_right[
        ['variant_read_count_propotion_per_replicate_left', 'variant_read_count_propotion_per_replicate_right']].apply(
        lambda row: row.min(), axis=1)


    distance_left_right = 1 - sum(variant_read_proportion_per_replicate_left_right['min_read_proportion'])
    return distance_left_right


def f12_filter_delete_renkonen(variant_read_count_df, renkonen_threshold):
    dfout = variant_read_count_df.copy()
    # dfout['filter_id'] = 12
    dfout['filter_delete'] = False
    #
    # group by on variant read count df  and aggregate by replicate_id to get all the replicate_id by biosample_id
    df2 = variant_read_count_df.groupby(['run_id', 'marker_id', 'biosample_id']).agg('replicate_id').apply(
        lambda x: list(set(x))).reset_index()
    df2['threshold_distance_number'] = df2['replicate_id'].apply(lambda x: (len(x) - 1) / 2)
    df2 = df2.loc[df2.threshold_distance_number != 0] # drop if threshold_distance_number == 0
    if df2.shape[0] == 0: # Only one replicate
        dfout['filter_delete'] = True
    else:
        df2['replicate_id_pairwise'] = df2.replicate_id.apply(lambda x: list(itertools.combinations(x, 2)))
        df2.drop('replicate_id', axis=1, inplace=True)
        df3 = pandas.DataFrame(
            data={'run_id': [], 'marker_id': [], 'biosample_id': [], 'left_replicate_id': [], 'right_replicate_id': [],
                  'renkonen_distance': []}, dtype='int')
        for row in df2.itertuples():
            run_id = row.run_id
            marker_id = row.marker_id
            biosample_id = row.biosample_id
            replicate_id_pairwise = row.replicate_id_pairwise
            for left_replicate_id, right_replicate_id in replicate_id_pairwise:
                df3 = pandas.concat(
                    [df3, pandas.DataFrame({'run_id': [run_id], 'marker_id': [marker_id], 'biosample_id': [biosample_id],
                                            'left_replicate_id': [left_replicate_id],
                                            'right_replicate_id': [right_replicate_id]})], axis=0, sort=True)
        # count the renkonen distance by pair of replicate for each biosample_id
        for row in df3.itertuples():
            run_id = row.run_id
            marker_id = row.marker_id
            biosample_id = row.biosample_id
            left_replicate_id = row.left_replicate_id
            right_replicate_id = row.right_replicate_id
            d = renkonen_distance(variant_read_count_df, run_id, marker_id, biosample_id, left_replicate_id,
                                  right_replicate_id)
            df3.loc[(df3.run_id == run_id) & (df3.marker_id == marker_id) & (df3.biosample_id == biosample_id)
                    & (df3.left_replicate_id == left_replicate_id) & (
                                df3.right_replicate_id == right_replicate_id), 'renkonen_distance'] = d
        # compare the renkonen distance to the renkonen_threshold
        df3['is_distance_gt_rthr'] = df3.renkonen_distance > renkonen_threshold
        # extract from the data frame df3 the combinaison of (replicate_left ,is_distance_gt_rthr) and (replicate_right ,is_distance_gt_rthr)
        df4 = pandas.DataFrame(
            data={'run_id': [], 'marker_id': [], 'biosample_id': [], 'replicate_id': [], 'is_distance_gt_rthr': []},
            dtype='int')
        df4 = df4.rename(columns={'replicate_id': 'left_replicate_id'})
        df4 = pandas.concat([df4, df3[['run_id', 'marker_id', 'biosample_id', 'left_replicate_id', 'is_distance_gt_rthr']]])
        df4 = df4.rename(columns={'left_replicate_id': 'right_replicate_id'})
        df4 = pandas.concat(
            [df4, df3[['run_id', 'marker_id', 'biosample_id', 'right_replicate_id', 'is_distance_gt_rthr']]], axis=0)
        df4 = df4.rename(columns={'right_replicate_id': 'replicate_id'})
        # group the data frame by 'run_id', 'marker_id', 'biosample_id', 'replicate_id' to count the sum  distance number for each replicate by biosample
        # merge with the df2 to get the threshold_distance_number
        df5 = df4.groupby(['run_id', 'marker_id', 'biosample_id', 'replicate_id']).sum().reset_index()
        df5 = df5.rename(columns={'is_distance_gt_rthr': 'distance_number'})
        dfout = dfout.merge(df2[['run_id', 'marker_id', 'biosample_id', 'threshold_distance_number']])
        dfout = dfout.merge(df5[['run_id', 'marker_id', 'biosample_id', 'replicate_id', 'distance_number']])
        #if  distance_number > threshold_distance_number do not pass the renkonen filter
        # df5['filter_delete'] = False
        dfout.loc[dfout.distance_number > dfout.threshold_distance_number, 'filter_delete'] = True
        #merge resulted data frame df5 with the variant_read_count_df
        # dfout = variant_read_count_df.merge(df5)
        dfout.drop(['distance_number', 'threshold_distance_number'], axis=1, inplace=True)
    #
    return dfout


