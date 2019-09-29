from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger

from wopmetabarcoding.utils.OptionManager import OptionManager
from sqlalchemy import select
import pandas


class OptimizeLFN(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.OptimizeLFN"
    }

    # Input file
    __input_file_positive_variants = "positive_variants"
    # Input table
    __input_table_run = "Run"
    __input_table_marker = "Marker"
    __input_table_biosample = "Biosample"
    __input_table_variant_read_count = "VariantReadCount"
    # Output file
    __output_file_optimize_lfn = "optimize_lfn"


    def specify_input_file(self):
        return[
            OptimizeLFN.__input_file_positive_variants,
        ]

    def specify_input_table(self):
        return [
            OptimizeLFN.__input_table_marker,
            OptimizeLFN.__input_table_run,
            OptimizeLFN.__input_table_biosample,
            OptimizeLFN.__input_table_variant_read_count,
        ]


    def specify_output_file(self):
        return [
            OptimizeLFN.__output_file_optimize_lfn,
        ]

    def specify_params(self):
        return {
            "log_verbosity": "int",
            "log_file": "str"
        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        OptionManager.instance()['log_verbosity'] = int(self.option("log_verbosity"))
        OptionManager.instance()['log_file'] = str(self.option("log_file"))

        ##########################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ##########################################################
        #
        # Input file path
        input_file_positive_variants = self.input_file(OptimizeLFN.__input_file_positive_variants)
        #
        # Input table models
        run_model = self.input_table(OptimizeLFN.__input_table_run)
        marker_model = self.input_table(OptimizeLFN.__input_table_marker)
        biosample_model = self.input_table(OptimizeLFN.__input_table_biosample)
        variant_read_count_model = self.input_table(OptimizeLFN.__input_table_variant_read_count)
        #
        # Output file path
        output_file_optimize_lfn = self.output_file(OptimizeLFN.__output_file_optimize_lfn)

        ##########################################################
        #
        # 1. Read input_file_positive_variants to get run_id, marker_id, biosample_id, for current analysis
        #
        ##########################################################
        positive_variant_df = pandas.read_csv(input_file_positive_variants, sep="\t", header=0,\
            names=['marker_name', 'run_name', 'biosample_name', 'variant_id', 'variant_sequence'], index_col=False)

        ##########################################################
        #
        # 2. Select run_id, marker_id, variant_id, biosample, replicate where variant, biosample, etc in positive_variants_df
        # 3. Get read_count: N_ijk
        #
        ##########################################################
        variant_read_count_list = []
        with engine.connect() as conn:
            for row in positive_variant_df.itertuples():
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
        optimized_lfn_df = pandas.DataFrame.from_records(variant_read_count_list,
            columns=['run_id', 'marker_id', 'variant_id', 'biosample_id', 'replicate_id', 'N_ijk'])
        optimized_lfn_df.drop_duplicates(inplace=True)

        ##########################################################
        #
        # 4. Compute ratio per_sum_biosample_replicate: N_ijk / N_jk
        #
        ##########################################################
        aggregate_df = optimized_lfn_df[['run_id', 'marker_id', 'biosample_id', 'replicate_id', 'N_ijk']].groupby(
            by=['run_id', 'marker_id', 'biosample_id', 'replicate_id']).sum().reset_index()
        aggregate_df = aggregate_df.rename(columns={'N_ijk': 'N_jk'})
        optimized_lfn_df = optimized_lfn_df.merge(aggregate_df, on=['run_id', 'marker_id', 'biosample_id', 'replicate_id'])
        optimized_lfn_df['lfn per_sum_biosample_replicate: N_ijk/N_jk'] = optimized_lfn_df['N_ijk'] / optimized_lfn_df['N_jk']

        ##########################################################
        #
        # 5. Compute ratio per_sum_variant: N_ijk / N_i
        #
        ##########################################################
        aggregate_df = optimized_lfn_df[['run_id', 'marker_id', 'variant_id', 'N_ijk']].groupby(
            by=['run_id', 'marker_id', 'variant_id']).sum().reset_index()
        aggregate_df = aggregate_df.rename(columns={'N_ijk': 'N_i'})
        optimized_lfn_df = optimized_lfn_df.merge(aggregate_df, on=['run_id', 'marker_id', 'variant_id'])
        optimized_lfn_df['lfn per_sum_variant: N_ijk/N_i'] = optimized_lfn_df['N_ijk'] / optimized_lfn_df['N_i']

        ##########################################################
        #
        # 6. Compute ratio per_sum_variant_replicate: N_ijk / N_ij
        #
        ##########################################################
        aggregate_df = optimized_lfn_df[['run_id', 'marker_id', 'variant_id', 'biosample_id', 'N_ijk']].groupby(
            by=['run_id', 'marker_id', 'variant_id', 'biosample_id']).sum().reset_index()
        aggregate_df = aggregate_df.rename(columns={'N_ijk': 'N_ij'})
        optimized_lfn_df = optimized_lfn_df.merge(aggregate_df, on=['run_id', 'marker_id', 'variant_id', 'biosample_id'])
        optimized_lfn_df['per_sum_variant_replicate: N_ijk/N_ij'] = optimized_lfn_df['N_ijk'] / optimized_lfn_df['N_ij']

        ##########################################################
        #
        # 7. Write TSV file
        #
        ##########################################################
        optimized_lfn_df.to_csv(output_file_optimize_lfn, header=True, sep='\t', float_format='%.10f', index=False)
