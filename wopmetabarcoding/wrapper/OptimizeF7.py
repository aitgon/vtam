import inspect
import os

from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger
from wopmetabarcoding.wrapper.FilterLFNutilities import FilterLFNRunner

# from wopmetabarcoding.wrapper.OptimizeLFNutilities import OptimizeLFNRunner
from sqlalchemy import select
import pandas

from wopmetabarcoding.utils.logger import logger


class OptimizeF7(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.OptimizeLFN"
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
    __output_file_optimize_lfn = "optimizef7"


    def specify_input_file(self):
        return[
            OptimizeF7.__input_file_positive_variants,
        ]

    def specify_input_table(self):
        return [
            OptimizeF7.__input_table_marker,
            OptimizeF7.__input_table_run,
            OptimizeF7.__input_table_biosample,
            OptimizeF7.__input_table_variant,
            OptimizeF7.__input_table_variant_read_count,
        ]


    def specify_output_file(self):
        return [
            OptimizeF7.__output_file_optimize_lfn,
        ]

    def specify_params(self):
        return {
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
        input_file_positive_variants = self.input_file(OptimizeF7.__input_file_positive_variants)
        #
        # Input table models
        run_model = self.input_table(OptimizeF7.__input_table_run)
        marker_model = self.input_table(OptimizeF7.__input_table_marker)
        biosample_model = self.input_table(OptimizeF7.__input_table_biosample)
        variant_model = self.input_table(OptimizeF7.__input_table_variant)

        variant_read_count_model = self.input_table(OptimizeF7.__input_table_variant_read_count)
        #
        # Output file path
        output_file_optimize_lfn = self.output_file(OptimizeF7.__output_file_optimize_lfn)


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
                   os.mknod(output_file_optimize_lfn)
                   exit()


        ##########################################################
        #
        # Select "keep" variants
        #
        ##########################################################
        logger.debug(
            "file: {}; line: {}; Extract some columns and 'keep' variants"
                .format(__file__, inspect.currentframe().f_lineno))
        variants_keep_df = variants_optimize_df[variants_optimize_df["action"].isin(["keep"])]
        variants_keep_df = variants_keep_df[['marker_name', 'run_name', 'biosample_name', 'biosample_type', 'variant_id', 'variant_sequence']]


        ##########################################################
        #
        # 2. Select run_id, marker_id, variant_id, biosample, replicate where variant, biosample, etc in positive_variants_df
        # 3. Get read_count: N_ijk
        #
        ##########################################################
        logger.debug(
            "file: {}; line: {}; Select run_id, marker_id, variant_id, biosample, replicate where variant, biosample, etc in positive_variants_df"
                .format(__file__, inspect.currentframe().f_lineno))
        logger.debug(
            "file: {}; line: {}; Get read_count: N_ijk"
                .format(__file__, inspect.currentframe().f_lineno))
        variant_read_count_keep_list = []
        with engine.connect() as conn:
            for row in variants_keep_df.itertuples():
                run_name = row.run_name
                marker_name = row.marker_name
                biosample_name = row.biosample_name
                variant_id = row.variant_id
                stmt_select = select([
                    run_model.__table__.c.id,
                    marker_model.__table__.c.id,
                    variant_read_count_model.__table__.c.variant_id,
                    biosample_model.__table__.c.id,
                    variant_read_count_model.__table__.c.replicate_id,
                    variant_read_count_model.__table__.c.read_count,])\
                    .where(run_model.__table__.c.name==run_name)\
                    .where(variant_read_count_model.__table__.c.run_id==run_model.__table__.c.id)\
                    .where(marker_model.__table__.c.name==marker_name)\
                    .where(variant_read_count_model.__table__.c.marker_id==marker_model.__table__.c.id)\
                    .where(biosample_model.__table__.c.name==biosample_name)\
                    .where(variant_read_count_model.__table__.c.biosample_id==biosample_model.__table__.c.id)\
                    .where(variant_read_count_model.__table__.c.variant_id==variant_id)\
                    .distinct()
                variant_read_count_keep_list = variant_read_count_keep_list + conn.execute(stmt_select).fetchall()
        variants_read_count_keep_df = pandas.DataFrame.from_records(variant_read_count_keep_list,
                                                       columns=['run_id', 'marker_id', 'variant_id', 'biosample_id',
                                                                  'replicate_id', 'N_ijk'])
        variants_read_count_keep_df.drop_duplicates(inplace=True)

        # ##########################################################
        # #
        # # N_ijk / N_ik
        # #
        # ##########################################################
        #
        # N_ik_df = variants_read_count_keep_df.groupby(by=['run_id', 'marker_id', 'variant_id', 'replicate_id']).sum().reset_index()
        # N_ik_df.drop(['biosample_id'], axis=1, inplace=True)
        # N_ik_df = N_ik_df.rename(columns={'N_ijk': 'N_ik'})
        # out_df = variants_read_count_keep_df.merge(N_ik_df, on=['run_id', 'marker_id', 'variant_id', 'replicate_id'])
        # out_df['N_ijk/N_ik'] = out_df.N_ijk / out_df.N_ik
        #
        # ##########################################################
        # #
        # # N_ijk / N_i
        # #
        # ###########
        # N_i_df = variants_read_count_keep_df.groupby(by=['run_id', 'marker_id', 'variant_id']).sum().reset_index()
        # N_i_df.drop(['biosample_id', 'replicate_id'], axis=1, inplace=True)
        # N_i_df = N_i_df.rename(columns={'N_ijk': 'N_i'})
        # out_df = out_df.merge(N_i_df, on=['run_id', 'marker_id', 'variant_id'])
        # out_df['N_ijk/N_i'] = out_df.N_ijk / out_df.N_i

        ##########################################################
        #
        # Select "delete" variants
        #
        ##########################################################

        variant_delete_df = variants_optimize_df[~variants_optimize_df["action"].isin(["keep", "tolerate"])]
        # variant_delete_df = variant_delete_df[['biosample_name', 'variant_id']]

        ##########################################################
        #
        # 2. Select run_id, marker_id, variant_id, biosample, replicate where variant, biosample, etc in positive_variants_df
        # 3. Get read_count: N_ijk
        #
        ##########################################################
        logger.debug(
            "file: {}; line: {}; Select run_id, marker_id, variant_id, biosample, replicate where variant, biosample, etc in positive_variants_df with action delete "
                .format(__file__, inspect.currentframe().f_lineno))
        logger.debug(
            "file: {}; line: {}; Get read_count: N_ijk"
                .format(__file__, inspect.currentframe().f_lineno))
        variant_read_count_delete_list = []

        with engine.connect() as conn:
            for row in variant_delete_df.itertuples():
                run_name = row.run_name
                marker_name = row.marker_name
                biosample_name = row.biosample_name
                # variant_id = row.variant_id
                stmt_select = select([
                    run_model.__table__.c.id,
                    marker_model.__table__.c.id,
                    variant_read_count_model.__table__.c.variant_id,
                    biosample_model.__table__.c.id,
                    variant_read_count_model.__table__.c.replicate_id,
                    variant_read_count_model.__table__.c.read_count, ]) \
                    .where(run_model.__table__.c.name == run_name) \
                    .where(variant_read_count_model.__table__.c.run_id == run_model.__table__.c.id) \
                    .where(marker_model.__table__.c.name == marker_name) \
                    .where(variant_read_count_model.__table__.c.marker_id == marker_model.__table__.c.id) \
                    .where(biosample_model.__table__.c.name == biosample_name) \
                    .where(variant_read_count_model.__table__.c.biosample_id == biosample_model.__table__.c.id) \
                    .distinct()
                variant_read_count_delete_list = variant_read_count_delete_list + conn.execute(stmt_select).fetchall()
        variants_read_count_delete_df = pandas.DataFrame.from_records(variant_read_count_delete_list,
                                                                    columns=['run_id', 'marker_id', 'variant_id',
                                                                             'biosample_id',
                                                                             'replicate_id', 'N_ijk'])
        variants_read_count_delete_df.drop_duplicates(inplace=True)

        ##########################################################
        #
        # Concatenate variants_read_count_keep_df and variants_read_count_delete_df
        #
        ##########################################################

        import pdb; pdb.set_trace()


        # ##########################################################
        # #
        # # N_ijk / N_ik for "delete" variants
        # #
        # ##########################################################
        #
        # N_ik_df = variants_read_count_delete_df.groupby(
        #     by=['run_id', 'marker_id', 'variant_id', 'replicate_id']).sum().reset_index()
        # N_ik_df.drop(['biosample_id'], axis=1, inplace=True)
        # N_ik_df = N_ik_df.rename(columns={'N_ijk': 'N_ik'})
        # out_delete_df = variants_read_count_delete_df.merge(N_ik_df, on=['run_id', 'marker_id', 'variant_id', 'replicate_id'])
        # out_delete_df['N_ijk/N_ik'] = out_delete_df.N_ijk / out_delete_df.N_ik
        #
        # ##########################################################
        # #
        # # N_ijk / N_i  for "delete" variants
        # #
        # ###########
        # N_i_df = variants_read_count_delete_df.groupby(by=['run_id', 'marker_id', 'variant_id']).sum().reset_index()
        # N_i_df.drop(['biosample_id', 'replicate_id'], axis=1, inplace=True)
        # N_i_df = N_i_df.rename(columns={'N_ijk': 'N_i'})
        # out_delete_df = out_delete_df.merge(N_i_df, on=['run_id', 'marker_id', 'variant_id'])
        # out_delete_df['N_ijk/N_i'] = out_delete_df.N_ijk / out_delete_df.N_i
        #
        #
        # ##########################################################
        # #
        # # 7. Write TSV file
        # #
        # ##########################################################
        # out_df['action'] = 'keep'
        #
        # # fo "delete" variants
        # out_delete_df['action'] = 'delete'
        # # concatenation
        # out_df.sort_values(by=['biosample_id', 'variant_id', 'N_ijk'], ascending=True, inplace=True)
        # out_delete_df.sort_values(by=['biosample_id', 'variant_id', 'N_ijk'], ascending=True, inplace=True)
        # out_df = pandas.concat([out_df, out_delete_df])
        #
        # # write to csv
        # out_df.to_csv(output_file_optimize_lfn, header=True, sep='\t', float_format='%.10f', index=False)
        #
        #
        #
        # ##########################################################
        # #
        # #   function of wrapper
        # #
        # ##########################################################
        #
        #  # column that should be on the variant_read_count_df  'id', 'run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id', 'read_count'
        #
        #  # run_id  marker_id  variant_id  biosample_id  replicate_id  N_ijk  N_ik  N_ijk/N_ik  N_i  N_ijk/N_i  action
        # variant_read_count_df = out_df[['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id', 'N_ijk']]
        # variant_read_count_df = variant_read_count_df.rename(columns={"N_ijk": "read_count"})
        #
        #
        # # Main loop to test parameters
        # # lfn_filter_runner = FilterLFNRunner(variant_read_count_df)
        # # import pdb;
        # # pdb.set_trace()
        # # lfn_filter_runner.f2_f4_lfn_delete_per_sum_variant(0.27325)
        # # lfn_filter_runner.f7_lfn_delete_absolute_read_count(63)
        # # lfn_filter_runner.delete_variant_df

