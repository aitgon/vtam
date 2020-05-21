import pandas
import sqlalchemy

from vtam.utils.SampleInformationFile import SampleInformationFile
from vtam.utils.KnownOccurrences import KnownOccurrences
from vtam.utils.VariantReadCountLikeDF import VariantReadCountLikeDF
from wopmars.models.ToolWrapper import ToolWrapper

from vtam.models.Run import Run
from vtam.models.Marker import Marker
from vtam.models.Biosample import Biosample
from vtam.models.Variant import Variant
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

        ##########################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ##########################################################

        # Input file output
        known_occurrences_tsv = self.input_file(
            OptimizeLFNbiosampleReplicate.__input_file_known_occurrences)
        fasta_info_tsv = self.input_file(
            OptimizeLFNbiosampleReplicate.__input_file_readinfo)

        # Output file output
        output_file_optimize_lfn = self.output_file(
            OptimizeLFNbiosampleReplicate.__output_file_optimize_lfn_biosample_replicate)

        ############################################################################################
        #
        # Get variant_read_count_df
        #
        ############################################################################################

        sample_info_tsv_obj = SampleInformationFile(tsv_path=fasta_info_tsv)
        variant_read_count_df = sample_info_tsv_obj.get_variant_read_count_df(
            VariantReadCount, engine=engine)

        ############################################################################################
        #
        # Get variant_read_count of "keep" variants in "mock" biosamples
        #
        ############################################################################################

        optimized_output_df = pandas.DataFrame()

        known_occurrences_df = KnownOccurrences(known_occurrences_tsv).to_identifier_df(engine)

        # Keep "mock" biosamples and "keep" variants
        known_occurrences_df = known_occurrences_df.loc[
            (known_occurrences_df.mock) & (known_occurrences_df.action == 'keep'), ]

        variant_read_count_df = variant_read_count_df.merge(known_occurrences_df, on=['run_id', 'marker_id', 'biosample_id',
                                                              'variant_id']).drop_duplicates()

        ############################################################################################
        #
        # Compute ratio per_sum_biosample_replicate: N_ijk / N_jk
        #
        ############################################################################################

        N_jk_df = VariantReadCountLikeDF(variant_read_count_df).get_N_jk_df()

        optimize_output_df = variant_read_count_df.merge(N_jk_df, on=[
            'run_id', 'marker_id', 'biosample_id', 'replicate'])
        optimize_output_df.rename(columns={'read_count': 'N_ijk'}, inplace=True)
        optimize_output_df['lfn_biosample_replicate: N_ijk/N_jk'] = optimize_output_df['N_ijk'] / \
            optimize_output_df['N_jk']

        ############################################################################################
        #
        # Sort and extract the lowest value of the ration with 4 digit
        #
        ############################################################################################

        optimize_output_df = optimize_output_df.sort_values(
            'lfn_biosample_replicate: N_ijk/N_jk', ascending=True)
        # Make round.inf with 4 decimals
        def round_down_4_decimals(x): return int(x * 10 ** 4) / 10 ** 4
        optimize_output_df['round_down'] = optimize_output_df['lfn_biosample_replicate: N_ijk/N_jk'].apply(
            round_down_4_decimals)

        ############################################################################################
        #
        # Add run, marker and biosample names
        #
        ############################################################################################

        run_df = pandas.read_sql(sqlalchemy.select([Run]), con=engine.connect())
        run_df.rename({'id': 'run_id', 'name': 'run', }, axis=1, inplace=True)
        optimize_output_df = optimize_output_df.merge(run_df, on='run_id')

        marker_df = pandas.read_sql(sqlalchemy.select([Marker]), con=engine.connect())
        marker_df.rename({'id': 'marker_id', 'name': 'marker', }, axis=1, inplace=True)
        optimize_output_df = optimize_output_df.merge(marker_df, on='marker_id')

        biosample_df = pandas.read_sql(sqlalchemy.select([Biosample]), con=engine.connect())
        biosample_df.rename({'id': 'biosample_id', 'name': 'biosample', }, axis=1, inplace=True)
        optimize_output_df = optimize_output_df.merge(biosample_df, on='biosample_id')

        ############################################################################################
        #
        # Compute ratios
        #
        ############################################################################################

        optimize_output_df['lfn_biosample_replicate: N_ijk/N_jk'] = optimize_output_df['N_ijk'] / optimize_output_df['N_jk']

        optimize_output_df = optimize_output_df.sort_values('lfn_biosample_replicate: N_ijk/N_jk', ascending=True)
        # Make round.inf with 4 decimals
        def round_down_4_decimals(x): return int(x * 10 ** 4) / 10 ** 4
        optimize_output_df['round_down'] = optimize_output_df['lfn_biosample_replicate: N_ijk/N_jk'].apply(round_down_4_decimals)

        ############################################################################################
        #
        # Rename columns, select columns and write to TSV file
        #
        ############################################################################################

        optimize_output_df = optimize_output_df[['run', 'marker', 'biosample', 'replicate', 'variant_id', 'N_ijk', 'N_jk', 'lfn_biosample_replicate: N_ijk/N_jk', 'round_down', 'variant_sequence']].drop_duplicates()
        optimize_output_df = optimize_output_df.rename({'variant_id': 'variant', 'variant_sequence': 'sequence'}, axis=1)

        #######################################################################
        #
        # 7. Write TSV file
        #
        #######################################################################

        optimize_output_df.sort_values(by=['run', 'marker', 'round_down', 'lfn_biosample_replicate: N_ijk/N_jk', 'biosample', 'replicate'], ascending=[True, True, True, True, True, True], inplace=True)
        optimize_output_df.to_csv(output_file_optimize_lfn, header=True, sep='\t',
            float_format='%.8f', index=False)

        # known_occurrences_grouped = known_occurrences_df.groupby(by=['run_id', 'marker_id'])
        # for run_marker_biosample_group in known_occurrences_grouped.groups:
        #     run_marker_row_df = known_occurrences_df.loc[
        #         known_occurrences_grouped.groups[run_marker_biosample_group],
        #         ['run_id', 'marker_id', 'biosample_id']]
        #     run_marker_row_df.drop_duplicates(inplace=True)
        #
        # #######################################################################
        # #
        # # Group and run_name this genetic_code by run_name/marker_name combination
        # # Loop by run_name/marker_name
        # #
        # #######################################################################
        #
        # # optimized_output_df = pandas.DataFrame()
        #
        # # known_occurrences_df = pandas.read_csv(
        # #     known_occurrences_tsv, sep="\t", header=0)
        # # known_occurrences_df.columns = known_occurrences_df.columns.str.lower()
        # # vknown_grouped = known_occurrences_df.groupby(by=['run', 'marker'])
        # for vknown_grouped_key in vknown_grouped.groups:
        #     run_marker_row_df = known_occurrences_df.loc[
        #         vknown_grouped.groups[vknown_grouped_key], :]
        #
        #     ###################################################################
        #     #
        #     # Read user known variant information and verify consistency with DB
        #     #
        #     ###################################################################
        #
        #     known_occurrences = KnownOccurrences(
        #         run_marker_row_df, fasta_info_tsv, engine)
        #     known_occurrences_ids_df = known_occurrences.known_occurrences_ids_df
        #
        #     ###################################################################
        #     #
        #     # Get run_id, marker_id, biosample_id rows marked as mock
        #     #
        #     ###################################################################
        #
        #     run_marker_biosample_mock_df = known_occurrences_ids_df.loc[
        #         known_occurrences_ids_df.biosample_type == 'mock', ['run_id', 'marker_id', 'biosample_id']]
        #     run_marker_biosample_mock_df.drop_duplicates(inplace=True)
        #
        #     ###################################################################
        #     #
        #     # Get variant read count variant_read_count_input_df
        #     #
        #     ###################################################################
        #
        #     # sample_info_tsv_obj = SampleInformationFile(
        #     #     tsv_path=fasta_info_tsv)
        #     # variant_read_count_df = sample_info_tsv_obj.get_variant_read_count_df(
        #     #     VariantReadCount, engine=engine)
        #
        #     variant_read_count_df_obj = VariantReadCountLikeDF(
        #         variant_read_count_df=variant_read_count_df)
        #     N_jk_df = variant_read_count_df_obj.get_N_jk_df()
        #
        #     ###################################################################
        #     #
        #     # Get delete variants, that are not keep in mock samples
        #     #
        #     ###################################################################
        #
        #     # These columns: run_id  marker_id  biosample_id  variant_id
        #     variant_keep_df = known_occurrences.get_run_marker_biosample_variant_df()
        #
        #     ###################################################################
        #     #
        #     # variant_read_count for keep variants only
        #     #
        #     ###################################################################
        #
        #     variant_read_count_keep_df = variant_read_count_df.merge(
        #         variant_keep_df, on=['run_id', 'marker_id', 'biosample_id', 'variant_id'])
        #
        #     ###################################################################
        #     #
        #     # 4. Compute ratio per_sum_biosample_replicate: N_ijk / N_jk
        #     #
        #     ###################################################################
        #
        #     optimize_output_df = variant_read_count_keep_df.merge(
        #         N_jk_df, on=['run_id', 'marker_id', 'biosample_id', 'replicate'])
        #     optimize_output_df.rename(
        #         columns={
        #             'read_count': 'N_ijk'},
        #         inplace=True)
        #     optimize_output_df['lfn_biosample_replicate: N_ijk/N_jk'] = optimize_output_df['N_ijk'] / \
        #         optimize_output_df['N_jk']
        #
        #     ###################################################################
        #     #
        #     # 5.Sort and extract the lowest value of the ration with 4 digit
        #     #
        #     ###################################################################
        #
        #     optimize_output_df = optimize_output_df.sort_values(
        #         'lfn_biosample_replicate: N_ijk/N_jk', ascending=True)
        #     # Make round.inf with 4 decimals
        #     def round_down_4_decimals(x): return int(x * 10 ** 4) / 10 ** 4
        #     optimize_output_df['round_down'] = optimize_output_df['lfn_biosample_replicate: N_ijk/N_jk'].apply(
        #         round_down_4_decimals)
        #
        #     ###################################################################
        #     #
        #     # Convert run_id, marker_id and biosample_id to their names
        #     #
        #     ###################################################################
        #
        #     with engine.connect() as conn:
        #         run_id_to_name = conn.execute(sqlalchemy.select(
        #             [Run.__table__.c.id, Run.__table__.c.name])).fetchall()
        #         marker_id_to_name = conn.execute(sqlalchemy.select(
        #             [Marker.__table__.c.id, Marker.__table__.c.name])).fetchall()
        #         biosample_id_to_name = conn.execute(sqlalchemy.select(
        #             [Biosample.__table__.c.id, Biosample.__table__.c.name])).fetchall()
        #
        #     run_id_to_name_df = pandas.DataFrame.from_records(
        #         data=run_id_to_name, columns=['run_id', 'run_name'])
        #     optimize_output_df = optimize_output_df.merge(
        #         run_id_to_name_df, on='run_id')
        #
        #     marker_id_to_name_df = pandas.DataFrame.from_records(
        #         data=marker_id_to_name, columns=['marker_id', 'marker_name'])
        #     optimize_output_df = optimize_output_df.merge(
        #         marker_id_to_name_df, on='marker_id')
        #
        #     biosample_id_to_name_df = pandas.DataFrame.from_records(
        #         data=biosample_id_to_name, columns=['biosample_id', 'biosample_name'])
        #     optimize_output_df = optimize_output_df.merge(
        #         biosample_id_to_name_df, on='biosample_id')
        #
        #     optimize_output_df = optimize_output_df[['run_name',
        #                                              'marker_name',
        #                                              'biosample_name',
        #                                              'replicate',
        #                                              'variant_id',
        #                                              'N_ijk',
        #                                              'N_jk',
        #                                              'lfn_biosample_replicate: N_ijk/N_jk',
        #                                              'round_down']]
        #
        #     optimized_output_df = pandas.concat(
        #         [optimized_output_df, optimize_output_df], axis=0)
        #
        # #######################################################################
        # #
        # # Add sequence
        # #
        # #######################################################################
        #
        # optimized_output_df.rename(
        #     {
        #         'marker_name': 'marker_name',
        #         'run_name': 'run_name',
        #         'biosample_name': 'biosample'},
        #     axis=1,
        #     inplace=True)
        # optimized_output_df['sequence'] = ''
        # with engine.connect() as conn:
        #     for variant_id in optimized_output_df.variant_id.unique():
        #         variant_id = int(variant_id)
        #         variant_sequence_row = conn.execute(sqlalchemy.select(
        #             [Variant.__table__.c.sequence]).where(Variant.__table__.c.id == variant_id)).first()
        #         if not (variant_sequence_row is None):
        #             variant_sequence = variant_sequence_row[0]
        #             optimized_output_df.loc[(
        #                 optimized_output_df.variant_id == variant_id).values, 'sequence'] = variant_sequence
        #
        # #######################################################################
        # #
        # # 7. Write TSV file
        # #
        # #######################################################################
        #
        # optimized_output_df.sort_values(
        #     by=[
        #         'run_name',
        #         'marker_name',
        #         'round_down',
        #         'lfn_biosample_replicate: N_ijk/N_jk',
        #         'biosample',
        #         'replicate'],
        #     ascending=[
        #         True,
        #         True,
        #         True,
        #         True,
        #         True,
        #         True],
        #     inplace=True)
        # optimized_output_df.to_csv(
        #     output_file_optimize_lfn,
        #     header=True,
        #     sep='\t',
        #     float_format='%.8f',
        #     index=False)
