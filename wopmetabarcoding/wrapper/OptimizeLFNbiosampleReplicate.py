import os
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger

# from wopmetabarcoding.wrapper.OptimizeLFNutilities import OptimizeLFNRunner
from sqlalchemy import select
import pandas

from wopmetabarcoding.utils.logger import logger


class OptimizeLFNbiosampleReplicate(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.OptimizeLFNbiosampleReplicate"
    }

    # Input file
    __input_file_variant_known = "variant_known"
    # Input table
    __input_table_run = "Run"
    __input_table_marker = "Marker"
    __input_table_biosample = "Biosample"
    __input_table_variant = "Variant"
    __input_table_variant_read_count = "VariantReadCount"
    # Output file
    __output_file_optimize_lfn = "optimize_lfn"


    def specify_input_file(self):
        return[
            OptimizeLFNbiosampleReplicate.__input_file_variant_known,
        ]

    def specify_input_table(self):
        return [
            OptimizeLFNbiosampleReplicate.__input_table_marker,
            OptimizeLFNbiosampleReplicate.__input_table_run,
            OptimizeLFNbiosampleReplicate.__input_table_biosample,
            OptimizeLFNbiosampleReplicate.__input_table_variant,
            OptimizeLFNbiosampleReplicate.__input_table_variant_read_count,
        ]


    def specify_output_file(self):
        return [
            OptimizeLFNbiosampleReplicate.__output_file_optimize_lfn,
        ]

    def specify_params(self):
        return {
        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind

        ##########################################################
        #
        # Wrapper inputs, outputs and parameters_numerical
        #
        ##########################################################
        #
        # Input file path
        input_file_variant_known = self.input_file(OptimizeLFNbiosampleReplicate.__input_file_variant_known)
        #
        # Input table models
        run_model = self.input_table(OptimizeLFNbiosampleReplicate.__input_table_run)
        marker_model = self.input_table(OptimizeLFNbiosampleReplicate.__input_table_marker)
        biosample_model = self.input_table(OptimizeLFNbiosampleReplicate.__input_table_biosample)
        variant_model = self.input_table(OptimizeLFNbiosampleReplicate.__input_table_variant)
        variant_read_count_model = self.input_table(OptimizeLFNbiosampleReplicate.__input_table_variant_read_count)
        #
        # Output file path
        output_file_optimize_lfn = self.output_file(OptimizeLFNbiosampleReplicate.__output_file_optimize_lfn)

        ##########################################################
        #
        # 1. Read variants_optimize to get run_id, marker_id, biosample_id, variant_id for current analysis
        #
        ##########################################################
        # positive_variant_df = pandas.read_csv(input_file_variant_known, sep="\t", header=0,\
        #     names=['marker_name', 'run_name', 'biosample_name', 'variant_id', 'variant_sequence'], index_col=False)
        variants_optimize_df = pandas.read_csv(input_file_variant_known, sep="\t", header=0, \
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
        # Extract some columns and "keep" variants
        #
        ##########################################################
        variants_keep_tolerate_df = variants_optimize_df[variants_optimize_df["action"].isin(["keep", "tolerate"])]
        variants_keep_tolerate_df = variants_keep_tolerate_df[['marker_name', 'run_name', 'biosample_name', 'variant_id', 'variant_sequence']]

        ##########################################################
        #
        # 2. Select run_id, marker_id, variant_id, biosample, replicate where variant, biosample, etc in variants_keep_df
        # 3. Get read_count: N_ijk
        #
        ##########################################################
        variant_read_count_list = []
        with engine.connect() as conn:
            for row in variants_keep_tolerate_df.itertuples():
                run_name = row.run_name
                marker_name = row.marker_name
                biosample_name = row.biosample_name
                variant_id = row.variant_id
                variant_sequence = row.variant_sequence
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
                variant_read_count_list = variant_read_count_list + conn.execute(stmt_select).fetchall()

        variant_read_count_df = pandas.DataFrame.from_records(variant_read_count_list,
            columns=['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id', 'N_ijk'])
        variant_read_count_df.drop_duplicates(inplace=True)


        ##########################################################
        #
        # 4. Compute ratio per_sum_biosample_replicate: N_ijk / N_jk
        #
        ##########################################################
        aggregate_df = variant_read_count_df[['run_id', 'marker_id', 'biosample_id', 'replicate_id', 'N_ijk']].groupby(
            by=['run_id', 'marker_id', 'biosample_id', 'replicate_id']).sum().reset_index()
        aggregate_df = aggregate_df.rename(columns={'N_ijk': 'N_jk'})
        variant_read_count_df = variant_read_count_df.merge(aggregate_df, on=['run_id', 'marker_id', 'biosample_id', 'replicate_id'])
        variant_read_count_df['lfn_biosample_replicate: N_ijk/N_jk'] = variant_read_count_df['N_ijk'] / variant_read_count_df['N_jk']

        ##########################################################
        #
        # 5.Sort and extract the lowest value of the ration with 4 digit
        #
        ##########################################################
        variant_read_count_df=variant_read_count_df.sort_values('lfn_biosample_replicate: N_ijk/N_jk', ascending=True)
        # Make round.inf with 4 decimals
        round_inf_4_decimals = lambda x: int(x * 10 ** 4) / 10 ** 4
        variant_read_count_df['round_inf'] = variant_read_count_df['lfn_biosample_replicate: N_ijk/N_jk'].apply(round_inf_4_decimals)

        ##########################################################
        #
        # 7. Write TSV file
        #
        ##########################################################
        variant_read_count_df.to_csv(output_file_optimize_lfn, header=True, sep='\t', float_format='%.10f', index=False)

