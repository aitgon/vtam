import pandas
import sqlalchemy

from vtam.models.Biosample import Biosample
from vtam.models.Marker import Marker
from vtam.models.Run import Run
from vtam.models.Variant import Variant
from vtam.utils.FilterPCRerrorRunner import FilterPCRerrorRunner
from vtam.utils.NameIdConverter import NameIdConverter


class OptimizePCRerrorRunner:

    def __init__(self, variant_read_count_df, known_occurrences_df):

        self.variant_read_count_df = variant_read_count_df
        self.known_occurrences_df = known_occurrences_df

    def get_optimize_df(self, engine):

        ############################################################################################
        #
        # Group and run_name this genetic_code by run_name/marker_name combination
        # Loop by run_name/marker_name
        #
        ############################################################################################

        optimize_df = pandas.DataFrame()

        variant_df = pandas.read_sql(sqlalchemy.select([
            Variant.__table__.c.id, Variant.__table__.c.sequence]), con=engine.connect(),
            index_col='id')

        known_occurrences_grouped = self.known_occurrences_df.groupby(
            by=['run_id', 'marker_id', 'biosample_id'])
        for run_marker_biosample_group in known_occurrences_grouped.groups:
            run_marker_biosample_df = self.known_occurrences_df.loc[
                known_occurrences_grouped.groups[run_marker_biosample_group],
                ['run_id', 'marker_id', 'biosample_id']]
            run_marker_biosample_df.drop_duplicates(inplace=True)

            variant_read_count_biosample_df = self.variant_read_count_df.merge(run_marker_biosample_df,
                                                                          on=['run_id', 'marker_id',
                                                                              'biosample_id'])

            ########################################################################################
            #
            # Get variant_expected_df, variant_unexpected_df
            #
            ########################################################################################

            variant_expected_df = (self.known_occurrences_df.merge(run_marker_biosample_df))[
                ['variant_id', 'variant_sequence']].drop_duplicates()
            variant_expected_df = pandas.DataFrame(
                {'sequence': variant_expected_df.variant_sequence.tolist()},
                index=variant_expected_df.variant_id.tolist())

            variant_unexpected_df = self.variant_read_count_df.merge(run_marker_biosample_df)[
                'variant_id'].drop_duplicates()
            variant_unexpected_df = variant_unexpected_df.loc[
                ~variant_unexpected_df.isin(variant_expected_df.index)]
            variant_unexpected_df = variant_df.loc[variant_unexpected_df.tolist(), :]

            ########################################################################################
            #
            # Run vsearch and get alignement variant_read_count_input_df
            #
            ########################################################################################

            filter_pcr_error_runner = FilterPCRerrorRunner(
                variant_expected_df=variant_expected_df, variant_unexpected_df=variant_unexpected_df,
                variant_read_count_df=variant_read_count_biosample_df)

            pcr_error_df = filter_pcr_error_runner.get_variant_unexpected_to_expected_ratio_df()

            optimize_df = pandas.concat(
                [optimize_df, pcr_error_df], axis=0)

        return optimize_df


    def to_tsv(self, optimize_path, engine):

        optimize_df = self.get_optimize_df(engine=engine)

        ############################################################################################
        #
        # Convert run_id, marker_id and biosample_id to their names
        #
        ############################################################################################

        optimize_df['run'] = NameIdConverter(optimize_df.run_id.tolist(),
                                                     engine).to_names(Run)
        optimize_df['marker'] = NameIdConverter(optimize_df.marker_id.tolist(),
                                                        engine).to_names(Marker)
        optimize_df['biosample'] = NameIdConverter(
            optimize_df.biosample_id.tolist(), engine).to_names(Biosample)
        optimize_df['sequence_expected'] = NameIdConverter(
            optimize_df.variant_id_expected.tolist(), engine).variant_id_to_sequence()
        optimize_df['sequence_unexpected'] = NameIdConverter(
            optimize_df.variant_id_unexpected.tolist(), engine).variant_id_to_sequence()

        optimize_df.drop(['run_id', 'marker_id', 'biosample_id'], axis=1, inplace=True)

        ############################################################################################
        #
        # Rename and order columns and write to TSV
        #
        ############################################################################################

        optimize_df.rename({'run_name': 'run_id',
                                    'marker_name': 'marker_id',
                                    'variant_id_expected': 'variant_expected',
                                    'variant_id_unexpected': 'variant_unexpected'}, axis=1,
                                   inplace=True)
        optimize_df = optimize_df[['run', 'marker', 'biosample', 'variant_expected',
                                                   'N_ij_expected', 'variant_unexpected',
                                                   'N_ij_unexpected',
                                                   'N_ij_unexpected_to_expected_ratio',
                                                   'sequence_expected',
                                                   'sequence_unexpected']]

        # lines should be ordered: by run_name, marker_name, ...
        optimize_df.sort_values(by=['run',
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

        optimize_df.to_csv(optimize_path, header=True, sep='\t', float_format='%.8f', index=False)
