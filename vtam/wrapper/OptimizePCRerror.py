import os
import pathlib

from vtam.models.VariantReadCount import VariantReadCount
from vtam.utils.KnownOccurrences import KnownOccurrences
from vtam.utils.OptimizePCRerrorRunner import OptimizePCRerrorRunner
from vtam.utils.PathManager import PathManager
from vtam.utils.SampleInformationFile import SampleInformationFile
from wopmars.models.ToolWrapper import ToolWrapper


class OptimizePCRerror(ToolWrapper):

    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.OptimizePCRerror"
    }

    # Input file
    __input_file_readinfo = "readinfo"
    __input_file_known_occurrences = "known_occurrences"
    # Input table
    __input_table_run = "Run"
    __input_table_marker = "Marker"
    __input_table_variant = "Variant"
    __input_table_biosample = "Biosample"
    __input_table_variant_read_count = "VariantReadCount"
    # Output file
    __output_file_optimize_pcr_error = "optimize_pcr_error"

    #

    def specify_input_file(self):
        return[
            OptimizePCRerror.__input_file_readinfo,
            OptimizePCRerror.__input_file_known_occurrences,
        ]

    def specify_input_table(self):
        return [
            OptimizePCRerror.__input_table_marker,
            OptimizePCRerror.__input_table_run,
            OptimizePCRerror.__input_table_variant,
            OptimizePCRerror.__input_table_biosample,
            OptimizePCRerror.__input_table_variant_read_count,
        ]

    def specify_output_file(self):
        return [
            OptimizePCRerror.__output_file_optimize_pcr_error,
        ]

    def specify_params(self):
        return {
        }

    def run(self):
        session = self.session
        engine = session._session().get_bind()

        this_temp_dir = os.path.join(
            PathManager.instance().get_tempdir(),
            os.path.basename(__file__))
        pathlib.Path(this_temp_dir).mkdir(exist_ok=True)

        ############################################################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ############################################################################################

        # Input file paths
        known_occurrences_tsv = self.input_file(OptimizePCRerror.__input_file_known_occurrences)
        fasta_info_tsv = self.input_file(OptimizePCRerror.__input_file_readinfo)
        #
        # Output file paths
        output_optimize_path = self.output_file(
            OptimizePCRerror.__output_file_optimize_pcr_error)

        ############################################################################################
        #
        # Get nijk_df, known_occurrences_df
        #
        ############################################################################################

        sample_info_tsv_obj = SampleInformationFile(tsv_path=fasta_info_tsv)
        variant_read_count_df = sample_info_tsv_obj.get_nijk_df(
            VariantReadCount, engine=engine)

        known_occurrences_df = KnownOccurrences(known_occurrences_tsv).to_identifier_df(engine)

        ############################################################################################
        #
        # Run optimizer and Write
        #
        ############################################################################################

        optimize_pcr_error_runner = OptimizePCRerrorRunner(
            variant_read_count_df=variant_read_count_df, known_occurrences_df=known_occurrences_df)
        optimize_pcr_error_runner.to_tsv(optimize_path=output_optimize_path, engine=engine)
