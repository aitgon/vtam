from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger

from wopmetabarcoding.wrapper.FilterLFNutilities import FilterLFNRunner
from sqlalchemy import select
import pandas


class OptimizePCRError(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.OptimizePCRError"
    }

    # Input file
    __input_file_positive_variants = "positive_variants"
    # Input table
    __input_table_run = "Run"
    __input_table_marker = "Marker"
    __input_table_biosample = "Biosample"
    __input_table_variant_read_count = "VariantReadCount"
    # Output file
    __output_file_optimize_lfn_per_biosample_replicate_threshold = "optimize_lfn_per_biosample_replicate_threshold"


    def specify_input_file(self):
        return[
            OptimizePCRError.__input_file_positive_variants,
        ]

    def specify_input_table(self):
        return [
            OptimizePCRError.__input_table_marker,
            OptimizePCRError.__input_table_run,
            OptimizePCRError.__input_table_biosample,
            OptimizePCRError.__input_table_variant_read_count,
        ]


    def specify_output_file(self):
        return [
            OptimizePCRError.__output_file_optimize_lfn_per_biosample_replicate_threshold,
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
        import pdb; pdb.set_trace()
        ################
        #
        # 1. Read positive_variants.tsv
        # 2. Select run_id, marker_id, variant_id, biosample, replicate where variant, biosample, etc in positive_variants_df
        # 3. Get read_count: N_ijk
        # 4. Compute ratio per_sum_biosample_replicate: N_ijk / N_jk
        # 5. Compute ratio per_sum_variant: N_ijk / N_i
        # 6. Compute ratio per_sum_variant_replicate: N_ijk / N_ij
        # 7. Output run_name, marker_name, variant, positive biosample j, replicate k. N_ijk, N_jk , N_ijk /N_jk, N_ijk / N_i, N_ijk / N_ij
        #
