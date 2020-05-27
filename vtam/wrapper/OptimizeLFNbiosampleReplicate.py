from vtam.utils.OptimizeLFNbiosampleReplicateRunner import OptimizeLFNbiosampleReplicateRunner
from vtam.utils.SampleInformationFile import SampleInformationFile
from vtam.utils.KnownOccurrences import KnownOccurrences
from wopmars.models.ToolWrapper import ToolWrapper
from vtam.models.VariantReadCount import VariantReadCount


class OptimizeLFNbiosampleReplicate(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.OptimizeLFNbiosampleReplicate"
    }

    # Input file
    __input_file_readinfo = "readinfo"
    __input_file_known_occurrences = "known_occurrences"
    # Input table
    __input_table_run = "Run"
    __input_table_marker = "Marker"
    __input_table_biosample = "Biosample"
    __input_table_variant = "Variant"
    __input_table_variant_read_count = "VariantReadCount"
    # Output file
    __output_file_optimize_lfn_biosample_replicate = "optimize_lfn_biosample_replicate"

    def specify_input_file(self):
        return[
            OptimizeLFNbiosampleReplicate.__input_file_readinfo,
            OptimizeLFNbiosampleReplicate.__input_file_known_occurrences,
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
            OptimizeLFNbiosampleReplicate.__output_file_optimize_lfn_biosample_replicate,
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
            OptimizeLFNbiosampleReplicate.__input_file_known_occurrences)
        fasta_info_tsv = self.input_file(
            OptimizeLFNbiosampleReplicate.__input_file_readinfo)

        # Output file output
        output_optimize_path = self.output_file(
            OptimizeLFNbiosampleReplicate.__output_file_optimize_lfn_biosample_replicate)

        ############################################################################################
        #
        # Get nijk_df and known_occurrences_df (keep)
        #
        ############################################################################################

        sample_info_tsv_obj = SampleInformationFile(tsv_path=fasta_info_tsv)
        variant_read_count_df = sample_info_tsv_obj.get_nijk_df(
            VariantReadCount, engine=engine)
        known_occurrences_df = KnownOccurrences(known_occurrences_tsv).to_identifier_df(engine)
        known_occurrences_df = known_occurrences_df.loc[
            (known_occurrences_df.mock == 1) & (known_occurrences_df.action == 'keep'), ]

        ############################################################################################
        #
        # Run optimizer and Write
        #
        ############################################################################################

        optimize_lfn_biosample_replicate_runner = OptimizeLFNbiosampleReplicateRunner(
            variant_read_count_df=variant_read_count_df, known_occurrences_df=known_occurrences_df)
        optimize_lfn_biosample_replicate_runner.to_tsv(optimize_path=output_optimize_path, engine=engine)


