from vtam.utils.RunnerOptimizeLFNsampleReplicate import RunnerOptimizeLFNsampleReplicate
from vtam.utils.FileSampleInformation import FileSampleInformation
from vtam.utils.FileKnownOccurrences import FileKnownOccurrences
from wopmars.models.ToolWrapper import ToolWrapper
from vtam.models.VariantReadCount import VariantReadCount


class OptimizeLFNsampleReplicate(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.OptimizeLFNsampleReplicate"
    }

    # Input file
    __input_file_sortedinfo = "sortedinfo"
    __input_file_known_occurrences = "known_occurrences"
    # Input table
    __input_table_run = "Run"
    __input_table_marker = "Marker"
    __input_table_sample = "Sample"
    __input_table_variant = "Variant"
    __input_table_variant_read_count = "VariantReadCount"
    # Output file
    __output_file_optimize_lfn_sample_replicate = "optimize_lfn_sample_replicate"

    def specify_input_file(self):
        return[
            OptimizeLFNsampleReplicate.__input_file_sortedinfo,
            OptimizeLFNsampleReplicate.__input_file_known_occurrences,
        ]

    def specify_input_table(self):
        return [
            OptimizeLFNsampleReplicate.__input_table_marker,
            OptimizeLFNsampleReplicate.__input_table_run,
            OptimizeLFNsampleReplicate.__input_table_sample,
            OptimizeLFNsampleReplicate.__input_table_variant,
            OptimizeLFNsampleReplicate.__input_table_variant_read_count,
        ]

    def specify_output_file(self):
        return [
            OptimizeLFNsampleReplicate.__output_file_optimize_lfn_sample_replicate,
        ]

    def specify_params(self):
        return {
        }

    def run(self):
        session = self.session
        engine = session._session().get_bind()

        ############################################################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ############################################################################################

        # Input file output
        known_occurrences_tsv = self.input_file(
            OptimizeLFNsampleReplicate.__input_file_known_occurrences)
        fasta_info_tsv = self.input_file(
            OptimizeLFNsampleReplicate.__input_file_sortedinfo)

        # Output file output
        output_optimize_path = self.output_file(
            OptimizeLFNsampleReplicate.__output_file_optimize_lfn_sample_replicate)

        ############################################################################################
        #
        # Get nijk_df and known_occurrences_df (keep)
        #
        ############################################################################################

        sample_info_tsv_obj = FileSampleInformation(tsv_path=fasta_info_tsv)
        variant_read_count_df = sample_info_tsv_obj.get_nijk_df(
            VariantReadCount, engine=engine)
        known_occurrences_df = FileKnownOccurrences(known_occurrences_tsv).to_identifier_df(engine)
        known_occurrences_df = known_occurrences_df.loc[
            (known_occurrences_df.mock == 1) & (known_occurrences_df.action == 'keep'), ]

        ############################################################################################
        #
        # Run optimizer and Write
        #
        ############################################################################################

        optimize_lfn_sample_replicate_runner = RunnerOptimizeLFNsampleReplicate(
            variant_read_count_df=variant_read_count_df, known_occurrences_df=known_occurrences_df)
        optimize_lfn_sample_replicate_runner.to_tsv(optimize_path=output_optimize_path, engine=engine)


