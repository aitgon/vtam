import os
import pathlib

import sqlalchemy
from vtam.utils.NameIdConverter import NameIdConverter
from wopmars.models.ToolWrapper import ToolWrapper

from vtam.utils.FilterPCRerrorRunner import FilterPCRerrorRunner
from vtam.utils.SampleInformationFile import SampleInformationFile
from vtam.utils.KnownOccurrences import KnownOccurrences
import pandas

from vtam.models.Run import Run
from vtam.models.Marker import Marker
from vtam.models.Biosample import Biosample
from vtam.models.Variant import Variant
from vtam.models.VariantReadCount import VariantReadCount

from vtam.utils.PathManager import PathManager


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
        output_file_optimize_pcr_error_path = self.output_file(
            OptimizePCRerror.__output_file_optimize_pcr_error)

        ############################################################################################
        #
        # Get variant_read_count_df and variant_df
        #
        ############################################################################################

        sample_info_tsv_obj = SampleInformationFile(tsv_path=fasta_info_tsv)
        variant_read_count_df = sample_info_tsv_obj.get_variant_read_count_df(
            VariantReadCount, engine=engine)

        variant_df = pandas.read_sql(sqlalchemy.select([
            Variant.__table__.c.id, Variant.__table__.c.sequence]), con=engine.connect(), index_col='id')

        ############################################################################################
        #
        # Group and run_name this genetic_code by run_name/marker_name combination
        # Loop by run_name/marker_name
        #
        ############################################################################################

        optimized_output_df = pandas.DataFrame()

        known_occurrences_df = KnownOccurrences(known_occurrences_tsv).to_identifier_df(engine)

        known_occurrences_df = known_occurrences_df.loc[known_occurrences_df.mock, ]

        known_occurrences_grouped = known_occurrences_df.groupby(by=['run_id', 'marker_id', 'biosample_id'])
        for run_marker_biosample_group in known_occurrences_grouped.groups:
            run_marker_biosample_df = known_occurrences_df.loc[
                known_occurrences_grouped.groups[run_marker_biosample_group],
                ['run_id', 'marker_id', 'biosample_id']]
            run_marker_biosample_df.drop_duplicates(inplace=True)

            variant_read_count_biosample_df = variant_read_count_df.merge(run_marker_biosample_df, on=['run_id', 'marker_id', 'biosample_id'])

            ############################################################################################
            #
            # Get variant_expected_df, variant_unexpected_df
            #
            ############################################################################################

            variant_expected_df = (known_occurrences_df.merge(run_marker_biosample_df))[
                ['variant_id', 'variant_sequence']].drop_duplicates()
            variant_expected_df = pandas.DataFrame({'sequence': variant_expected_df.variant_sequence.tolist()},
                             index=variant_expected_df.variant_id.tolist())

            variant_unexpected_df = variant_read_count_df.merge(run_marker_biosample_df)['variant_id'].drop_duplicates()
            variant_unexpected_df = variant_unexpected_df.loc[~variant_unexpected_df.isin(variant_expected_df.index)]
            variant_unexpected_df = variant_df.loc[variant_unexpected_df.tolist(), :]

            ########################################################################################
            #
            # Run vsearch and get alignement variant_read_count_input_df
            #
            ########################################################################################

            filter_pcr_error_runner = FilterPCRerrorRunner(
                variant_expected_df=variant_expected_df,
                variant_unexpected_df=variant_unexpected_df,
                variant_read_count_df=variant_read_count_biosample_df)

            pcr_error_df = filter_pcr_error_runner.get_variant_unexpected_to_expected_ratio_df()

            optimized_output_df = pandas.concat(
                [optimized_output_df, pcr_error_df], axis=0)

        ############################################################################################
        #
        # Convert run_id, marker_id and biosample_id to their names
        #
        ############################################################################################

        optimized_output_df['run'] = NameIdConverter(optimized_output_df.run_id.tolist(), engine).to_names(Run)
        optimized_output_df['marker'] = NameIdConverter(optimized_output_df.marker_id.tolist(), engine).to_names(Marker)
        optimized_output_df['biosample'] = NameIdConverter(optimized_output_df.biosample_id.tolist(), engine).to_names(Biosample)
        optimized_output_df['sequence_expected'] = NameIdConverter(optimized_output_df.variant_id_expected.tolist(), engine).variant_id_to_sequence()
        optimized_output_df['sequence_unexpected'] = NameIdConverter(optimized_output_df.variant_id_unexpected.tolist(), engine).variant_id_to_sequence()

        optimized_output_df.drop(['run_id', 'marker_id', 'biosample_id'], axis=1, inplace=True)

        ############################################################################################
        #
        # Rename and order columns and write to TSV
        #
        ############################################################################################

        optimized_output_df.rename({'run_name': 'run_id',
                                    'marker_name': 'marker_id',
                                    'variant_id_expected': 'variant_expected',
                                    'variant_id_unexpected': 'variant_unexpected'}, axis=1, inplace=True)
        optimized_output_df = optimized_output_df[['run', 'marker', 'biosample', 'variant_expected',
                                                   'N_ij_expected', 'variant_unexpected',
                                                   'N_ij_unexpected',
                                                   'N_ij_unexpected_to_expected_ratio',
                                                   'sequence_expected',
                                                   'sequence_unexpected']]

        ############################################################################################
        #
        # Write Optimize PCRError to TSV file
        #
        ############################################################################################

        # lines should be ordered: by run_name, marker_name, ...
        optimized_output_df.sort_values(by=['run',
                                           'marker',
                                           'N_ij_unexpected_to_expected_ratio',
                                           'variant_expected',
                                           'variant_unexpected'],
                                       ascending=[True,
                                                  True,
                                                  False,
                                                  True,
                                                  True],
                                       inplace=True)
        optimized_output_df.to_csv(
            output_file_optimize_pcr_error_path,
            header=True,
            sep='\t',
            float_format='%.10f',
            index=False)
