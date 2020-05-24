import pandas
import sqlalchemy

from vtam.models.Biosample import Biosample
from vtam.models.Marker import Marker
from vtam.models.Run import Run
from vtam.utils.VariantReadCountLikeDF import VariantReadCountLikeDF


class OptimizeLFNbiosampleReplicateRunner:

    def __init__(self, variant_read_count_df, known_occurrences_df):

        self.variant_read_count_df = variant_read_count_df
        self.known_occurrences_df = known_occurrences_df

    def get_optimize_df(self, engine):

        ############################################################################################
        #
        # OptimizeLFNbiosampleReplicateRunner
        #
        ############################################################################################

        variant_read_count_df = self.variant_read_count_df.merge(self.known_occurrences_df, on=['run_id', 'marker_id', 'biosample_id',
                                                              'variant_id']).drop_duplicates()

        ############################################################################################
        #
        # Compute ratio per_sum_biosample_replicate: N_ijk / N_jk
        #
        ############################################################################################

        N_jk_df = VariantReadCountLikeDF(variant_read_count_df).get_N_jk_df()

        optimize_df = variant_read_count_df.merge(N_jk_df, on=[
            'run_id', 'marker_id', 'biosample_id', 'replicate'])
        optimize_df.rename(columns={'read_count': 'N_ijk'}, inplace=True)
        optimize_df['lfn_biosample_replicate: N_ijk/N_jk'] = optimize_df['N_ijk'] / \
            optimize_df['N_jk']

        ############################################################################################
        #
        # Sort and extract the lowest value of the ration with 4 digit
        #
        ############################################################################################

        optimize_df = optimize_df.sort_values(
            'lfn_biosample_replicate: N_ijk/N_jk', ascending=True)
        # Make round.inf with 4 decimals
        def round_down_4_decimals(x): return int(x * 10 ** 4) / 10 ** 4
        optimize_df['round_down'] = optimize_df['lfn_biosample_replicate: N_ijk/N_jk'].apply(
            round_down_4_decimals)

        ############################################################################################
        #
        # Add run, marker and biosample names
        #
        ############################################################################################

        run_df = pandas.read_sql(sqlalchemy.select([Run]), con=engine.connect())
        run_df.rename({'id': 'run_id', 'name': 'run', }, axis=1, inplace=True)
        optimize_df = optimize_df.merge(run_df, on='run_id')

        marker_df = pandas.read_sql(sqlalchemy.select([Marker]), con=engine.connect())
        marker_df.rename({'id': 'marker_id', 'name': 'marker', }, axis=1, inplace=True)
        optimize_df = optimize_df.merge(marker_df, on='marker_id')

        biosample_df = pandas.read_sql(sqlalchemy.select([Biosample]), con=engine.connect())
        biosample_df.rename({'id': 'biosample_id', 'name': 'biosample', }, axis=1, inplace=True)
        optimize_df = optimize_df.merge(biosample_df, on='biosample_id')

        ############################################################################################
        #
        # Compute ratios
        #
        ############################################################################################

        optimize_df['lfn_biosample_replicate: N_ijk/N_jk'] = optimize_df['N_ijk'] / optimize_df['N_jk']

        optimize_df = optimize_df.sort_values('lfn_biosample_replicate: N_ijk/N_jk', ascending=True)
        # Make round.inf with 4 decimals
        def round_down_4_decimals(x): return int(x * 10 ** 4) / 10 ** 4
        optimize_df['round_down'] = optimize_df['lfn_biosample_replicate: N_ijk/N_jk'].apply(round_down_4_decimals)

        ############################################################################################
        #
        # Rename columns, select columns and write to TSV file
        #
        ############################################################################################

        optimize_df = optimize_df[['run', 'marker', 'biosample', 'replicate', 'variant_id', 'N_ijk', 'N_jk', 'lfn_biosample_replicate: N_ijk/N_jk', 'round_down', 'variant_sequence']].drop_duplicates()
        optimize_df = optimize_df.rename({'variant_id': 'variant', 'variant_sequence': 'sequence'}, axis=1)

        #######################################################################
        #
        # 7. Write TSV file
        #
        #######################################################################

        optimize_df.sort_values(by=['run', 'marker', 'round_down', 'lfn_biosample_replicate: N_ijk/N_jk', 'biosample', 'replicate'], ascending=[True, True, True, True, True, True], inplace=True)

        return optimize_df

    def to_tsv(self, optimize_path, engine):

        optimize_df = self.get_optimize_df(engine)

        optimize_df.to_csv(optimize_path, header=True, sep='\t',
            float_format='%.8f', index=False)
