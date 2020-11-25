import numpy
from vtam.utils.RunnerOptimizeLFNreadCountAndVariantRunMarker import \
    RunnerOptimizeLFNreadCountAndVariantRunMarker

from vtam.utils.constants import lfn_ni_njk_cutoff_global_max, lfn_nijk_cutoff_global_max, \
    lfn_nijk_cutoff_lst_size

from vtam.models.Marker import Marker
from vtam.models.Run import Run
from vtam.models.VariantReadCount import VariantReadCount
from vtam.utils.FileKnownOccurrences import FileKnownOccurrences
from vtam.utils.NameIdConverter import NameIdConverter
from vtam.utils.RunnerOptimizeLFNreadCountAndVariant import \
    RunnerOptimizeLFNreadCountAndVariant
from vtam.utils.FileSampleInformation import FileSampleInformation
from wopmars.models.ToolWrapper import ToolWrapper


class OptimizeLFNreadCountAndLFNvariant(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.OptimizeLFNreadCountAndLFNvariant"
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
    __output_file_optimize_lfn_read_count_and_lfn_variant = "optimize_lfn_read_count_and_lfn_variant"
    __output_file_optimize_lfn_variant_specific = "optimize_lfn_variant_specific"

    def specify_input_file(self):
        return [
            OptimizeLFNreadCountAndLFNvariant.__input_file_sortedinfo,
            OptimizeLFNreadCountAndLFNvariant.__input_file_known_occurrences,
            "params",
        ]

    def specify_input_table(self):
        return [
            OptimizeLFNreadCountAndLFNvariant.__input_table_marker,
            OptimizeLFNreadCountAndLFNvariant.__input_table_run,
            OptimizeLFNreadCountAndLFNvariant.__input_table_sample,
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
            "lfn_variant_cutoff": "float",
            "lfn_variant_replicate_cutoff": "float",
            "lfn_sample_replicate_cutoff": "required|float",
            "lfn_read_count_cutoff": "required|float",
            "min_replicate_number": "required|int",
        }

    def run(self):
        """
        Algorithm (Updated Oct 13, 2019)

        1. Read file with known variants (Mock/tolerate, delete and real)
        2. Control if user variants and sequence are consistent in the database
        3. Get variant_read_count of this run_name-marker_name-sample-replicate experiment
        5. Compute maximal lfn_nijk_cutoff that keeps all 'keep' variants with the 'run_lfn_read_count_and_lfn_variant' algorithm
        6. Compute maximal lfn_variant_cutoff that keeps all 'keep' variants with the 'run_lfn_read_count_and_lfn_variant' algorithm (See below)
        7. Loop between default and lfn_nijk_cutoff and run_lfn_read_count_and_lfn_variant parameters
            7.1 Compute number of keep variants. Should be always maximal.
            7.2 Compute number of delete variants Should decrease.
        8. Compute variant(-replicate) specific cutoff for delete variants
            8.1 For each variant i (Or variant-replicate i-k ),
                get N_ijk_max and use it to computer variant specific cutoff

        Description of the 'run_lfn_read_count_and_lfn_variant' algorithm

        1. Remove if does not pass these filter
            1.1 Filter lfn_variant (Or lfn_variant_replicate)
            1.2 Filter lfn_sample_replicate
            1.3 Filter absolute read count
        2. Filter if not min replicate number

        """
        session = self.session
        engine = session._session().get_bind()

        ############################################################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ############################################################################################

        # Input file output
        known_occurrences_tsv = self.input_file(OptimizeLFNreadCountAndLFNvariant.__input_file_known_occurrences)
        fasta_info_tsv = self.input_file(OptimizeLFNreadCountAndLFNvariant.__input_file_sortedinfo)

        # Output file output
        output_file_optimize_lfn_tsv = self.output_file(
            OptimizeLFNreadCountAndLFNvariant.__output_file_optimize_lfn_read_count_and_lfn_variant)
        output_file_lfn_variant_specific_cutoff_tsv = self.output_file(
            OptimizeLFNreadCountAndLFNvariant.__output_file_optimize_lfn_variant_specific)

        # Options
        lfn_ni_cutoff = self.option("lfn_variant_cutoff")
        lfn_nik_cutoff = self.option("lfn_variant_replicate_cutoff")
        min_replicate_number = self.option("min_replicate_number")
        lfn_njk_cutoff = self.option("lfn_sample_replicate_cutoff")
        lfn_nijk_cutoff = int(self.option("lfn_read_count_cutoff"))

        filter_kwargs = {"lfn_ni_cutoff": lfn_ni_cutoff,
                         "lfn_nik_cutoff": lfn_nik_cutoff,
                         "lfn_njk_cutoff": lfn_njk_cutoff,
                         "lfn_nijk_cutoff": lfn_nijk_cutoff,
                         'min_replicate_number': min_replicate_number,
                      }

        ############################################################################################
        #
        # Get nijk_df and known_occurrences_df (keep)
        #
        ############################################################################################

        sample_info_tsv_obj = FileSampleInformation(tsv_path=fasta_info_tsv)
        nijk_df = sample_info_tsv_obj.get_nijk_df(
            VariantReadCount, engine=engine)

        known_occurrences_df = FileKnownOccurrences(known_occurrences_tsv).to_identifier_df(engine)

        ############################################################################################
        #
        # Create cutoff values lists
        #
        ############################################################################################

        # # lfn_nijk_cutoff_list = range(lfn_nijk_cutoff, lfn_nijk_cutoff_global_max + 1, round(int((lfn_nijk_cutoff_global_max - lfn_nijk_cutoff + 1)/10), -1))
        # lfn_nijk_cutoff_list = range(lfn_nijk_cutoff, lfn_nijk_cutoff_global_max + 1, round(int((lfn_nijk_cutoff_global_max - lfn_nijk_cutoff + 1)/10), -1))
        # lfn_nijk_cutoff_list = RunnerOptimizeLFNreadCountAndVariantRunMarker.get_lfn_nijk_cutoff_lst(start=lfn_nijk_cutoff, stop=lfn_nijk_cutoff_global_max, nb_points=10)
        # lfn_nijk_cutoff_list = RunnerOptimizeLFNreadCountAndVariantRunMarker.get_lfn_nijk_cutoff_lst(start=lfn_nijk_cutoff, stop=lfn_nijk_cutoff_global_max, nb_points=10)
        # if lfn_nik_cutoff is None:  # lfn_variant optimization
        #     lfn_ni_nik_cutoff_list = [round(x, 3) for x in numpy.arange(lfn_ni_cutoff, lfn_ni_njk_cutoff_global_max + 0.001, (lfn_ni_njk_cutoff_global_max - lfn_ni_cutoff + 0.001)/10)]
        # else:  # lfn_variant_replicate optimization
        #     lfn_ni_nik_cutoff_list = [round(x, 3) for x in numpy.arange(lfn_ni_cutoff, lfn_ni_njk_cutoff_global_max + 0.001, (lfn_ni_njk_cutoff_global_max - lfn_ni_cutoff + 0.001)/10)]

        ############################################################################################
        #
        # Group and run_name this genetic_code by run_name/marker_name combination
        # Loop by run_name/marker_name
        #
        ############################################################################################

        optim_lfn_readcount_variant_runner = RunnerOptimizeLFNreadCountAndVariant(
            nijk_df=nijk_df, known_occurrences_df=known_occurrences_df)
        out_optimize_df, out_optimize2_df = optim_lfn_readcount_variant_runner.get_optimize_df(
            lfn_ni_cutoff=lfn_ni_cutoff, lfn_nik_cutoff=lfn_nik_cutoff, lfn_njk_cutoff=lfn_njk_cutoff,
            lfn_nijk_cutoff=lfn_nijk_cutoff, min_replicate_number=min_replicate_number)

        ############################################################################################
        #
        # out_optimize_df: Format and write
        #
        ############################################################################################

        out_optimize_df.marker_id = NameIdConverter(out_optimize_df.marker_id, engine=engine).to_names(Marker)
        out_optimize_df.run_id = NameIdConverter(out_optimize_df.run_id, engine=engine).to_names(Run)
        out_optimize_df.rename({'run_id': 'run', 'marker_id': 'marker'}, axis=1, inplace=True)
        out_optimize_df.to_csv(output_file_optimize_lfn_tsv, header=True, sep='\t', index=False)

        ############################################################################################
        #
        # out_optimize_df: Format and write
        #
        ############################################################################################

        out_optimize2_df.marker_id = NameIdConverter(out_optimize2_df.marker_id, engine=engine).to_names(Marker)
        out_optimize2_df.run_id = NameIdConverter(out_optimize2_df.run_id, engine=engine).to_names(Run)
        out_optimize2_df['action'] = 'delete'
        out_optimize2_df['sequence'] = NameIdConverter(out_optimize2_df.variant_id, engine=engine).variant_id_to_sequence()
        out_optimize2_df.rename({'run_id': 'run', 'marker_id': 'marker', 'variant_id': 'variant', 'read_count': 'read_count_max'}, axis=1, inplace=True)

        if self.option("lfn_variant_replicate_cutoff") is None:
            out_optimize2_df = out_optimize2_df[['run', 'marker', 'variant', 'action', 'read_count_max', 'N_i', 'lfn_variant_cutoff', 'sequence']]
        else:
            out_optimize2_df = out_optimize2_df[['run', 'marker', 'variant', 'replicate', 'action', 'read_count_max', 'N_ik', 'lfn_variant_replicate_cutoff', 'sequence']]

        out_optimize2_df.to_csv(
            output_file_lfn_variant_specific_cutoff_tsv, header=True, sep='\t', index=False)
