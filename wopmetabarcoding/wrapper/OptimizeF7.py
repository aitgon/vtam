import inspect
import os

from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger

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
        # 1. Read variants_optimize to get run_id, marker_id, biosample_id, variant_id for current analysis
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
        # Extract some columns and not in "keep or tolerate " variants
        #
        ##########################################################
        variant_delete_df = variants_optimize_df[~variants_optimize_df["action"].isin(["keep", "tolerate"])]
        variant_delete_df = variant_delete_df[['biosample_name', 'variant_id']]


        ##########################################################
        #
        # Extract some columns and "keep" variants
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
        #  3. Get read_count: N_ijk
        #
        ##########################################################
        logger.debug(
            "file: {}; line: {}; Select run_id, marker_id, variant_id, biosample, replicate where variant, biosample, etc in positive_variants_df"
                .format(__file__, inspect.currentframe().f_lineno))
        logger.debug(
            "file: {}; line: {}; Get read_count: N_ijk"
                .format(__file__, inspect.currentframe().f_lineno))
        variant_read_count_list = []
        with engine.connect() as conn:
            for row in variant_delete_df.itertuples():
                biosample_name = row.biosample_name
                variant_id = row.variant_id
                stmt_select = select([
                    run_model.__table__.c.id,
                    marker_model.__table__.c.id,
                    variant_read_count_model.__table__.c.variant_id,
                    biosample_model.__table__.c.id,
                    variant_read_count_model.__table__.c.replicate_id,
                    variant_read_count_model.__table__.c.read_count, ]) \
                    .where(biosample_model.__table__.c.name != biosample_name) \
                    .where(variant_read_count_model.__table__.c.variant_id != variant_id) \
                    .distinct()
                variant_read_count_list = variant_read_count_list + conn.execute(stmt_select).fetchall()
        variants_read_count_df = pandas.DataFrame.from_records(variant_read_count_list,
                                                       columns=['run_id', 'marker_id', 'variant_id', 'biosample_id',
                                                                  'replicate_id', 'N_ijk'])
        variants_read_count_df.drop_duplicates(inplace=True)


        ##########################################################
        #
        # 4 - Create data frame containing the variant_id, lfn_read_count_threshold, sample_type (mock1,mock2 ...)
        #
        ##########################################################
        logger.debug(
            "file: {}; line: {}; Create data frame containing the variant_id, lfn_read_count_threshold, sample_type (mock1,mock2 ...)"
                .format(__file__, inspect.currentframe().f_lineno))
        lfn_read_count_threshod_dic = {}
        lfn_read_count_threshod_list = []

        biosample_type_list = variants_keep_df.biosample_type.tolist()
        for biosample_type in biosample_type_list:
            df = variants_keep_df[variants_keep_df["biosample_type"] == biosample_type]
            df = df[['marker_name', 'run_name', 'biosample_name','biosample_type','variant_id', 'variant_sequence']]
            variant_id_list = df.variant_id.tolist()
            for variant_id in variant_id_list:
                df_by_biosample_type = variants_read_count_df[variants_read_count_df["variant_id"] == variant_id]
                rdc_list =df_by_biosample_type["N_ijk"].tolist()
                rdc_max_min  = [max(rdc_list), min(rdc_list)]
                if  rdc_list[1] not in rdc_max_min: lfn_read_count_threshold = rdc_list[1]
                if rdc_list[2] not in rdc_max_min: lfn_read_count_threshold = rdc_list[2]
                if rdc_list[3] not in rdc_max_min: lfn_read_count_threshold = rdc_list[3]
                #
                read_count_number= df_by_biosample_type["N_ijk"].sum()
                lfn_read_count_threshod_dic = {}
                # write in data frame the variant_id with his lfn_read_count_threshold and sample_type
                lfn_read_count_threshod_dic["variant_id"] = variant_id
                import pdb; pdb.set_trace()
                lfn_read_count_threshod_dic["lfn_read_count_threshold"] = lfn_read_count_threshold
                lfn_read_count_threshod_dic["biosample_type"] = biosample_type
                lfn_read_count_threshod_dic["read_count_number"] = read_count_number
            lfn_read_count_threshod_list.append(lfn_read_count_threshod_dic)
        # data frame containing the ['variant_id','lfn_read_count_threshold','biosample_type','read_count_number']
        lfn_read_count_threshod_df = pandas.DataFrame(lfn_read_count_threshod_list)

        #choose the lowest value of the lfn_read_count_threshold from the data frame
        lfn_read_count_threshod_df = lfn_read_count_threshod_df.sort_values('lfn_read_count_threshold', ascending=True)
        lfn_read_count_threshold = lfn_read_count_threshod_df['lfn_read_count_threshold'].head(1).values[0]

        # for the lfn_per_variant_threshold Do the ration  lfn_read_count_threshol / number_read_count and choose alsoo the lowest one

        lfn_read_count_threshod_df["lfn_per_variant_threshold"] = lfn_read_count_threshod_df["lfn_read_count_threshold"]/ lfn_read_count_threshod_df["read_count_number"]
        # sort by lowest value lfn_per_variant_threshold
        lfn_read_count_threshod_df = lfn_read_count_threshod_df.sort_values('lfn_per_variant_threshold', ascending=True)
        #
        lfn_per_variant_threshold = lfn_read_count_threshod_df['lfn_per_variant_threshold'].head(1).values[0]


        ##########################################################
        #
        # 2. Select run_id, marker_id, variant_id, biosample, replicate where variant, biosample, etc in positive_variants_df
        #  3. Get read_count: N_ijk
        #
        ##########################################################
        logger.debug(
            "file: {}; line: {}; Select run_id, marker_id, variant_id, biosample, replicate where variant, biosample, etc in positive_variants_df"
                .format(__file__, inspect.currentframe().f_lineno))
        logger.debug(
            "file: {}; line: {}; Get read_count: N_ijk"
                .format(__file__, inspect.currentframe().f_lineno))
        variant_read_count_list = []
        with engine.connect() as conn:
            for row in variant_delete_df.itertuples():
                biosample_name = row.biosample_name
                variant_id = row.variant_id
                stmt_select = select([
                    run_model.__table__.c.id,
                    marker_model.__table__.c.id,
                    variant_read_count_model.__table__.c.variant_id,
                    biosample_model.__table__.c.id,
                    variant_read_count_model.__table__.c.replicate_id,
                    variant_read_count_model.__table__.c.read_count, ]) \
                    .where(biosample_model.__table__.c.name != biosample_name) \
                    .where(variant_read_count_model.__table__.c.variant_id != variant_id) \
                    .distinct()
                variant_read_count_list = variant_read_count_list + conn.execute(stmt_select).fetchall()
        variants_read_count_df = pandas.DataFrame.from_records(variant_read_count_list,
                                                       columns=['run_id', 'marker_id', 'variant_id', 'biosample_id',
                                                                  'replicate_id', 'N_ijk'])
        variants_read_count_df.drop_duplicates(inplace=True)




        ##########################################################
        #
        # 5- Run filter for a series of (lfn_read_count_threshold -lfn_per_variant_threshold) combinations
        #
        ##########################################################
        import pdb;
        pdb.set_trace()
        # for lfn_per_variant_threshold in itertools.product([1, 2], [5, 6], ['eleven', 'f']):
        #     for lfn_read_count_threshold in itertools.product([1, 2], [5, 6], ['eleven', 'f']):





        ##########################################################
        #
        # 6 -determine a specific lfn_per_variant_threshold specified to the delete variant that steel there despite using the optimal combinaison in 5)-
        #
        ##########################################################



        ##########################################################
        #
        # 7. Write TSV file
        #
        ##########################################################
        # optimized_lfn_df.to_csv(output_file_optimize_lfn, header=True, sep='\t', float_format='%.10f', index=False)
