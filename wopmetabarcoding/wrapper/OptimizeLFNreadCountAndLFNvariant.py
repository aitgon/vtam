import inspect
import os

import math
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger
from wopmetabarcoding.wrapper.FilterLFNutilities import FilterLFNRunner

from sqlalchemy import select
import pandas

from wopmetabarcoding.utils.logger import logger
from wopmetabarcoding.wrapper.FilterMinReplicateNumber import f9_delete_min_replicate_number


class OptimizeLFNreadCountAndLFNvariant(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.OptimizeLFNreadCountAndLFNvariant"
    }

    # Input file
    __input_file_positive_variants = "positive_variants"
    # Input table
    __input_table_run = "Run"
    __input_table_marker = "Marker"
    __input_table_biosample = "Biosample"
    __input_table_variant = "Variant"
    __input_table_variant_read_count = "VariantReadCount"
    # Output file
    __output_file_optimize_lfn_read_count_and_lfn_variant = "optimize_lfn_read_count_and_lfn_variant"
    __output_file_optimize_lfn_variant_specific = "optimize_lfn_variant_specific"


    def specify_input_file(self):
        return[
            OptimizeLFNreadCountAndLFNvariant.__input_file_positive_variants,
        ]

    def specify_input_table(self):
        return [
            OptimizeLFNreadCountAndLFNvariant.__input_table_marker,
            OptimizeLFNreadCountAndLFNvariant.__input_table_run,
            OptimizeLFNreadCountAndLFNvariant.__input_table_biosample,
            OptimizeLFNreadCountAndLFNvariant.__input_table_variant,
            OptimizeLFNreadCountAndLFNvariant.__input_table_variant_read_count,
        ]


    def specify_output_file(self):
        return [
            OptimizeLFNreadCountAndLFNvariant.__output_file_optimize_lfn_read_count_and_lfn_variant,
            OptimizeLFNreadCountAndLFNvariant.__output_file_optimize_lfn_variant_specific,
        ]

    def specify_params(self):
        return {
            "lfn_variant_threshold": "float",
            "lfn_biosample_replicate_threshold": "float",
            "lfn_read_count_threshold": "float",
            "min_replicate_number": "int",
        }

    def run(self, f7_lfn_delete_absolute_read_count=None):
        session = self.session()
        engine = session._WopMarsSession__session.bind

        ##########################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ##########################################################
        #
        # Input file path
        input_file_positive_variants = self.input_file(OptimizeLFNreadCountAndLFNvariant.__input_file_positive_variants)
        #
        # Input table models
        run_model = self.input_table(OptimizeLFNreadCountAndLFNvariant.__input_table_run)
        marker_model = self.input_table(OptimizeLFNreadCountAndLFNvariant.__input_table_marker)
        biosample_model = self.input_table(OptimizeLFNreadCountAndLFNvariant.__input_table_biosample)
        variant_model = self.input_table(OptimizeLFNreadCountAndLFNvariant.__input_table_variant)

        variant_read_count_model = self.input_table(OptimizeLFNreadCountAndLFNvariant.__input_table_variant_read_count)
        #
        # Output file path
        output_file_optimize_lfn_tsv = self.output_file(
            OptimizeLFNreadCountAndLFNvariant.__output_file_optimize_lfn_read_count_and_lfn_variant)
        output_file_lfn_variant_specific_threshold_tsv = self.output_file(
            OptimizeLFNreadCountAndLFNvariant.__output_file_optimize_lfn_variant_specific)
        #
        # Options
        min_replicate_number = self.option("min_replicate_number")
        lfn_biosample_replicate_threshold = self.option("lfn_biosample_replicate_threshold")
        lfn_read_count_threshold = self.option("lfn_read_count_threshold")
        lfn_variant_threshold = self.option("lfn_variant_threshold")

        ##########################################################
        #
        # Read "variants_optimize" to get run_id, marker_id, biosample_id, variant_id for current analysis
        #
        ##########################################################
        # positive_variant_df = pandas.read_csv(input_file_positive_variants, sep="\t", header=0,\
        #     names=['marker_name', 'run_name', 'biosample_name', 'variant_id', 'variant_sequence'], index_col=False)

        variants_optimize_df = pandas.read_csv(input_file_positive_variants, sep="\t", header=0, \
                                              names=['marker_name', 'run_name', 'biosample_name', 'biosample_type',
                                                     'variant_id', 'action', 'variant_sequence', 'note'], index_col=False)

        ########################
        #
        # Control if user variants and sequence are consistent in the database
        #
        ########################

        variant_control_df = variants_optimize_df[['variant_id', 'variant_sequence']].drop_duplicates()
        variant_control_df = variant_control_df.loc[~variant_control_df.variant_id.isnull()]
        with engine.connect() as conn:
            for row in variant_control_df.itertuples():
                variant_id = row.variant_id
                variant_sequence = row.variant_sequence
                stmt_select = select([variant_model.__table__.c.id, variant_model.__table__.c.sequence])\
                    .where(variant_model.__table__.c.id == variant_id)\
                    .where(variant_model.__table__.c.sequence == variant_sequence)
                if conn.execute(stmt_select).first() is None:
                   logger.error("Variant {} and its sequence are not coherent with the VTAM database".format(variant_id))
                   os.mknod(output_file_optimize_lfn_tsv)
                   exit()


        ##########################################################
        #
        # Select "keep" variants
        #
        ##########################################################
        # logger.debug(
        #     "file: {}; line: {}; Extract some columns and 'keep' variants"
        #         .format(__file__, inspect.currentframe().f_lineno))
        # variants_keep_df = variants_optimize_df[variants_optimize_df["action"].isin(["keep"])]
        # variants_keep_df = variants_keep_df[['marker_name', 'run_name', 'biosample_name', 'biosample_type', 'variant_id', 'variant_sequence']]

        ##########################################################
        #
        # Keep: Select run_id, marker_id, variant_id, biosample, replicate, N_ijk where variant, biosample, etc in positive_variants_df
        #
        ##########################################################
        logger.debug(
            "file: {}; line: {}; Select run_id, marker_id, variant_id, biosample, replicate, N_ijk where variant, biosample, etc in positive_variants_df"
                .format(__file__, inspect.currentframe().f_lineno))
        variant_read_count_list = []
        with engine.connect() as conn:
            for row in variants_optimize_df.itertuples():
                run_name = row.run_name
                marker_name = row.marker_name
                biosample_name = row.biosample_name
                biosample_type = row.biosample_type
                action = row.action
                stmt_select = select([
                    run_model.__table__.c.id,
                    marker_model.__table__.c.id,
                    variant_read_count_model.__table__.c.variant_id,
                    biosample_model.__table__.c.id,
                    variant_read_count_model.__table__.c.replicate_id,
                    variant_read_count_model.__table__.c.read_count])\
                    .where(run_model.__table__.c.name==run_name)\
                    .where(variant_read_count_model.__table__.c.run_id==run_model.__table__.c.id)\
                    .where(marker_model.__table__.c.name==marker_name)\
                    .where(variant_read_count_model.__table__.c.marker_id==marker_model.__table__.c.id)\
                    .where(biosample_model.__table__.c.name==biosample_name)\
                    .where(variant_read_count_model.__table__.c.biosample_id==biosample_model.__table__.c.id)

                if not math.isnan(row.variant_id): # variant id defined
                    variant_id = row.variant_id
                    stmt_select = stmt_select.where(variant_read_count_model.__table__.c.variant_id == variant_id)
                stmt_select = stmt_select.distinct()

                stmt_select_fetchall = conn.execute(stmt_select).fetchall()
                # append action
                variant_read_count_list = variant_read_count_list\
                                          + [list(r) + [biosample_type, action] for r in stmt_select_fetchall]

        variant_read_count_df = pandas.DataFrame.from_records(variant_read_count_list,
                                                       columns=['run_id', 'marker_id', 'variant_id', 'biosample_id',
                                                                  'replicate_id', 'read_count', 'biosample_type', 'action'])
        variant_read_count_df.drop_duplicates(inplace=True)

        ##########################################################
        #
        # Get keep, delete-negative and delete-real
        #
        ##########################################################

        variant_keep_df = variant_read_count_df.loc[variant_read_count_df.action == 'keep']
        # variant_keep_df.drop(['replicate_id', 'read_count', 'biosample_type', 'action'], axis=1, inplace=True)
        variant_keep_df.drop_duplicates(inplace=True)

        variant_delete_negative_df = variant_read_count_df.loc[(variant_read_count_df.action == 'delete') &
                                                     (variant_read_count_df.biosample_type == 'negative')]
        # variant_delete_negative_df.drop(['replicate_id', 'read_count', 'biosample_type', 'action'], axis=1, inplace=True)
        variant_delete_negative_df.drop_duplicates(inplace=True)


        variant_delete_real_df = variant_read_count_df.loc[(variant_read_count_df.action == 'delete') &
                                                     (variant_read_count_df.biosample_type == 'real')]
        # variant_delete_real_df.drop(['replicate_id', 'read_count', 'biosample_type', 'action'], axis=1, inplace=True)
        variant_delete_real_df.drop_duplicates(inplace=True)

        ##############
        #
        # Set maximal value of lfn_read_count_threshold
        #
        ##############
        #
        #
        count_keep = 0
        count_keep_max = 0
        #
        variant_read_count_df = variant_read_count_df[
            ['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id', 'read_count']]
        lfn_read_count_threshold_previous = 10
        #  loop over lfn_read_count_threshold
        for lfn_read_count_threshold in list(range(lfn_read_count_threshold_previous, 1001, 10)):

            variant_read_count_remained_df, count_keep = lfn_read_count_and_lfn_variant(variant_read_count_df, variant_keep_df, lfn_variant_threshold,
                                           lfn_biosample_replicate_threshold,
                                           lfn_read_count_threshold, min_replicate_number)

            if count_keep > count_keep_max:
                count_keep_max = count_keep
            elif count_keep < count_keep_max:
                break # stop when count_keep starts to decrease
            lfn_read_count_threshold_previous = lfn_read_count_threshold

        lfn_read_count_threshold_max = lfn_read_count_threshold_previous # upper border of lfn_read_count_threshold

        ##############
        #
        # Set maximal value of lfn_variant_threshold
        #
        ##############
        #

        count_keep = 0
        count_keep_max = 0
        #
        variant_read_count_df = variant_read_count_df[
            ['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id', 'read_count']]
        lfn_variant_threshold_previous = 0.001 # is divided by 1000
        # loop over lfn_variant_threshold: 0.001, 0.002, ...
        for lfn_variant_threshold in [i/1000 for i in range(int(lfn_variant_threshold_previous*1000), 101, 1)]:

            variant_read_count_remained_df, count_keep = lfn_read_count_and_lfn_variant(variant_read_count_df, variant_keep_df, lfn_variant_threshold,
                                           lfn_biosample_replicate_threshold,
                                           lfn_read_count_threshold, min_replicate_number)

            if count_keep > count_keep_max:
                count_keep_max = count_keep
            elif count_keep < count_keep_max:
                break # stop when count_keep starts to decrease
            lfn_variant_threshold_previous = lfn_variant_threshold

        lfn_variant_threshold_max = lfn_variant_threshold_previous # upper border of lfn_read_count_threshold

        ##############
        #
        # Optimize together two parameters to previous define borders
        #
        ##############
        #
        out_lfn_variant_list = []

        variant_read_count_df = variant_read_count_df[
            ['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id', 'read_count']]
        # lfn_read_count_threshold_default = lfn_read_count_threshold
        # lfn_variant_threshold_default = lfn_variant_threshold
        # while count_keep >= count_keep_max:
        for lfn_read_count_threshold in list(range(10, lfn_read_count_threshold_max, 10)): # loop over lfn_read_count_threshold
            # loop over lfn_variant_threshold: 0.001, 0.002, ...
            for lfn_variant_threshold in [i/1000 for i in range(1, int(lfn_variant_threshold_max*1000), 1)]:

                variant_read_count_remained_df, count_keep = lfn_read_count_and_lfn_variant(variant_read_count_df, variant_keep_df, lfn_variant_threshold,
                                               lfn_biosample_replicate_threshold,
                                               lfn_read_count_threshold, min_replicate_number)

                ##########################################################
                #
                # Count delete
                #
                ##########################################################

                variant_read_count_remained_delete_negative_df = variant_read_count_remained_df.merge(
                    variant_delete_negative_df, on=['run_id', 'marker_id', 'biosample_id']).drop_duplicates()
                variant_read_count_remained_delete_real_df = variant_read_count_remained_df.merge(
                    variant_delete_real_df, on=['run_id', 'marker_id', 'variant_id', 'biosample_id']).drop_duplicates()
                count_delete = variant_read_count_remained_delete_negative_df.shape[0] \
                               + variant_read_count_remained_delete_real_df.shape[0]

                ##########################################################
                #
                # Store results
                #
                ##########################################################
                if count_keep >= count_keep_max: # Store results if count_keep maximal
                    out_lfn_variant_row_dic = {"lfn_variant_threshold": lfn_variant_threshold,
                               "lfn_read_count_threshold": lfn_read_count_threshold,
                               "variant_nb_keep": count_keep, "variant_nb_delete": count_delete}
                    out_lfn_variant_list.append(out_lfn_variant_row_dic)
                # print(len(out_lfn_variant_list))

                #
                # if count_keep < count_keep_max:
                #     break
                if count_keep > count_keep_max:
                    count_keep_max = count_keep

        ##########################################################
        #
        # Write TSV file
        #
        ##########################################################

        out_lfn_variant_df = pandas.DataFrame(out_lfn_variant_list, columns=['variant_nb_keep', 'variant_nb_delete',
                                                                             'lfn_read_count_threshold',
                                                                             'lfn_variant_threshold'])

        out_lfn_variant_df.sort_values(by=["variant_nb_keep", "variant_nb_delete", "lfn_read_count_threshold",
                                           "lfn_variant_threshold"], ascending=[False, True, False, False], inplace=True)

        # out_lfn_variant_df.to_csv(output_file_optimize_lfn_tsv, header=True, sep='\t', float_format='%.10f', index=False)
        out_lfn_variant_df.to_csv(output_file_optimize_lfn_tsv, header=True, sep='\t', index=False)

        ##########################################################
        #
        # LFN variant specific threshold
        #
        ##########################################################

        variant_delete_df = pandas.concat([variant_delete_negative_df, variant_delete_real_df], axis=0)
        lfn_variant_specific_threshold_df = variant_delete_df.copy()
        lfn_variant_specific_threshold_df.drop(['action', 'biosample_type'], axis=1, inplace=True)
        N_i_df = lfn_variant_specific_threshold_df.groupby(by=['run_id', 'marker_id', 'variant_id'])\
            .sum().reset_index()
        N_i_df = N_i_df[['run_id', 'marker_id', 'variant_id', 'read_count']]
        N_i_df = N_i_df.rename(columns={'read_count': 'N_i'})
        N_i_df = N_i_df.drop_duplicates()
        lfn_variant_specific_threshold_df = lfn_variant_specific_threshold_df.merge(N_i_df, on=['run_id', 'marker_id', 'variant_id'])
        lfn_variant_specific_threshold_df['lfn_variant_threshold'] = lfn_variant_specific_threshold_df.read_count / lfn_variant_specific_threshold_df.N_i
        lfn_variant_specific_threshold_df.sort_values(by='lfn_variant_threshold', ascending=False, inplace=True)
        lfn_variant_specific_threshold_df.drop_duplicates('variant_id', keep='first', inplace=True)
        lfn_variant_specific_threshold_df = (lfn_variant_specific_threshold_df[
            ['run_id', 'marker_id', 'variant_id', 'read_count', 'N_i', 'lfn_variant_threshold']]).drop_duplicates()
        lfn_variant_specific_threshold_df = lfn_variant_specific_threshold_df.merge(variant_delete_df[
                                           ['run_id', 'marker_id', 'variant_id', 'biosample_type',
                                            'action']].drop_duplicates(), on=['run_id', 'marker_id', 'variant_id'])
        # lfn_variant_specific_threshold_df.to_csv(output_file_lfn_variant_specific_threshold_tsv, header=True, sep='\t', float_format='%.10f', index=False)
        lfn_variant_specific_threshold_df.to_csv(output_file_lfn_variant_specific_threshold_tsv, header=True, sep='\t', index=False)


def lfn_read_count_and_lfn_variant(variant_read_count_df, variant_keep_df, lfn_variant_threshold, lfn_biosample_replicate_threshold,
                                   lfn_read_count_threshold, min_replicate_number):
    lfn_filter_runner = FilterLFNRunner(variant_read_count_df)

    ###################
    #
    # Filter lfn_variant
    #
    ####################

    lfn_filter_runner.f2_f4_lfn_delete_variant(lfn_variant_threshold)

    ###################
    #
    # Filter lfn_biosample_replicate
    #
    ####################

    lfn_filter_runner.f6_lfn_delete_biosample_replicate(lfn_biosample_replicate_threshold)

    ###################
    #
    # Filter absolute read count
    #
    ####################

    lfn_filter_runner.f7_lfn_delete_absolute_read_count(lfn_read_count_threshold)

    ###################
    #
    # f8_lfn_delete_do_not_pass_all_filters
    #
    ####################

    lfn_filter_runner.f8_lfn_delete_do_not_pass_all_filters()

    variant_read_count_remained_df = lfn_filter_runner.delete_variant_df

    variant_read_count_remained_df = variant_read_count_remained_df.loc[
        (variant_read_count_remained_df.filter_id == 8) &
        (variant_read_count_remained_df.filter_delete == 0)]

    ##########################################################
    #
    # f9_delete_min_replicate_number
    #
    ##########################################################

    variant_read_count_remained_df = f9_delete_min_replicate_number(variant_read_count_remained_df, min_replicate_number)
    variant_read_count_remained_df = variant_read_count_remained_df.loc[
        (variant_read_count_remained_df.filter_delete == 0)]
    variant_read_count_remained_df.drop('filter_delete', axis=1, inplace=True)

    ##########################################################
    #
    # Count keep
    #
    ##########################################################

    variant_read_count_remained_df = variant_read_count_remained_df[['run_id', 'marker_id', 'variant_id', 'biosample_id']]
    variant_read_count_remained_df.drop_duplicates(inplace=True)
    variant_read_count_remained_keep_df = variant_read_count_remained_df.merge(variant_keep_df,
                                                                               on=['run_id', 'marker_id', 'variant_id',
                                                                                   'biosample_id'])
    variant_read_count_remained_keep_df = variant_read_count_remained_keep_df[
        ['run_id', 'marker_id', 'variant_id', 'biosample_id']].drop_duplicates()
    count_keep = variant_read_count_remained_keep_df.shape[0]

    # Delete object

    del (lfn_filter_runner)

    return variant_read_count_remained_df, count_keep