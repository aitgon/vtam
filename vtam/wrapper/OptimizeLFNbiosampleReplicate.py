import pandas
import sqlalchemy
from wopmars.models.ToolWrapper import ToolWrapper

from vtam.utils.SampleInformationUtils import FastaInformationTSV
from vtam.utils.VariantKnown import VariantKnown
from vtam.utils.VariantReadCountDF import VariantReadCountDF


class OptimizeLFNbiosampleReplicate(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.OptimizeLFNbiosampleReplicate"
    }

    # Input file
    __input_file_readinfo = "readinfo"
    __input_file_variant_known = "variant_known"
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
            OptimizeLFNbiosampleReplicate.__input_file_variant_known,
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
        variant_known_tsv = self.input_file(OptimizeLFNbiosampleReplicate.__input_file_variant_known)
        fasta_info_tsv_path = self.input_file(OptimizeLFNbiosampleReplicate.__input_file_readinfo)
        #
        # Input table models
        run_model = self.input_table(OptimizeLFNbiosampleReplicate.__input_table_run)
        marker_model = self.input_table(OptimizeLFNbiosampleReplicate.__input_table_marker)
        biosample_model = self.input_table(OptimizeLFNbiosampleReplicate.__input_table_biosample)
        variant_model = self.input_table(OptimizeLFNbiosampleReplicate.__input_table_variant)
        variant_read_count_model = self.input_table(OptimizeLFNbiosampleReplicate.__input_table_variant_read_count)
        #
        # Output file output
        output_file_optimize_lfn = self.output_file(OptimizeLFNbiosampleReplicate.__output_file_optimize_lfn_biosample_replicate)

        ################################################################################################################
        #
        # Group and run this code by run/marker combination
        # Loop by run/marker
        #
        ################################################################################################################

        final_optimize_output_df = pandas.DataFrame()

        variant_known_df = pandas.read_csv(variant_known_tsv, sep="\t", header=0)
        variant_known_df.columns = variant_known_df.columns.str.lower()
        vknown_grouped = variant_known_df.groupby(by=['run', 'marker'])
        for vknown_grouped_key in vknown_grouped.groups:
            variant_known_i_df = variant_known_df.loc[vknown_grouped.groups[vknown_grouped_key], :]

            ################################################################################################################
            #
            # Read user known variant information and verify consistency with DB
            #
            ########################### #####################################################################################

            variant_known = VariantKnown(variant_known_i_df, fasta_info_tsv_path, engine)
            variant_known_ids_df = variant_known.variant_known_ids_df

            ################################################################################################################
            #
            # Get run_id, marker_id, biosample_id rows marked as mock
            #
            ################################################################################################################

            run_marker_biosample_mock_df = variant_known_ids_df.loc[
                variant_known_ids_df.biosample_type == 'mock', ['run_id', 'marker_id', 'biosample_id']]
            run_marker_biosample_mock_df.drop_duplicates(inplace=True)

            ################################################################################################################
            #
            # Get variant read count df
            #
            ################################################################################################################
            fasta_info_tsv_obj = FastaInformationTSV(fasta_info_tsv=fasta_info_tsv_path, engine=engine)
            variant_read_count_df = fasta_info_tsv_obj.get_variant_read_count_df(variant_read_count_model)

            variant_read_count_df_obj = VariantReadCountDF(variant_read_count_df=variant_read_count_df)
            N_jk_df = variant_read_count_df_obj.get_N_jk_df()

            ################################################################################################
            #
            # Get delete variants, that are not keep in mock samples
            #
            ################################################################################################

            # These columns: run_id  marker_id  biosample_id  variant_id
            variant_keep_df = variant_known.get_keep_run_marker_biosample_variant_df()


            ################################################################################################################
            #
            # variant_read_count for keep variants only
            #
            ################################################################################################################

            variant_read_count_keep_df = variant_read_count_df.merge(variant_keep_df, on=['run_id', 'marker_id',
                                                                                                        'biosample_id',
                                                                                                        'variant_id'])

            ##########################################################
            #
            # 4. Compute ratio per_sum_biosample_replicate: N_ijk / N_jk
            #
            ##########################################################

            optimize_output_df = variant_read_count_keep_df.merge(N_jk_df, on=['run_id', 'marker_id', 'biosample_id', 'replicate'])
            optimize_output_df.rename(columns={'read_count': 'N_ijk'}, inplace=True)
            optimize_output_df['lfn_biosample_replicate: N_ijk/N_jk'] = optimize_output_df['N_ijk'] / optimize_output_df['N_jk']

            ##########################################################
            #
            # 5.Sort and extract the lowest value of the ration with 4 digit
            #
            ##########################################################

            optimize_output_df=optimize_output_df.sort_values('lfn_biosample_replicate: N_ijk/N_jk', ascending=True)
            # Make round.inf with 4 decimals
            round_down_4_decimals = lambda x: int(x * 10 ** 4) / 10 ** 4
            optimize_output_df['round_down'] = optimize_output_df['lfn_biosample_replicate: N_ijk/N_jk'].apply(round_down_4_decimals)

            ##########################################################
            #
            # Convert run_id, marker_id and biosample_id to their names
            #
            ##########################################################

            with engine.connect() as conn:
                run_id_to_name = conn.execute(
                    sqlalchemy.select([run_model.__table__.c.id, run_model.__table__.c.name])).fetchall()
                marker_id_to_name = conn.execute(
                    sqlalchemy.select([marker_model.__table__.c.id, marker_model.__table__.c.name])).fetchall()
                biosample_id_to_name = conn.execute(
                    sqlalchemy.select([biosample_model.__table__.c.id, biosample_model.__table__.c.name])).fetchall()

            run_id_to_name_df = pandas.DataFrame.from_records(data=run_id_to_name, columns=['run_id', 'run_name'])
            optimize_output_df = optimize_output_df.merge(run_id_to_name_df, on='run_id')

            marker_id_to_name_df = pandas.DataFrame.from_records(data=marker_id_to_name, columns=['marker_id', 'marker_name'])
            optimize_output_df = optimize_output_df.merge(marker_id_to_name_df, on='marker_id')

            biosample_id_to_name_df = pandas.DataFrame.from_records(data=biosample_id_to_name, columns=['biosample_id', 'biosample_name'])
            optimize_output_df = optimize_output_df.merge(biosample_id_to_name_df, on='biosample_id')

            optimize_output_df = optimize_output_df[['run_name', 'marker_name', 'biosample_name', 'replicate', 'variant_id',
           'N_ijk', 'N_jk', 'lfn_biosample_replicate: N_ijk/N_jk', 'round_down']]

            final_optimize_output_df = pandas.concat([final_optimize_output_df, optimize_output_df], axis=0)

        ################################################################################################################
        #
        # Add sequence
        #
        ################################################################################################################

        final_optimize_output_df.rename({'marker_name': 'marker', 'run_name': 'run', 'biosample_name': 'biosample'}, axis=1, inplace=True)
        final_optimize_output_df['sequence'] = ''
        with engine.connect() as conn:
            for variant_id in final_optimize_output_df.variant_id.unique():
                variant_id = int(variant_id)
                variant_sequence_row = conn.execute(sqlalchemy.select([variant_model.__table__.c.sequence]).where(variant_model.__table__.c.id == variant_id)).first()
                if not (variant_sequence_row is None):
                    variant_sequence = variant_sequence_row[0]
                    final_optimize_output_df.loc[(final_optimize_output_df.variant_id == variant_id).values, 'sequence'] = variant_sequence

        ##########################################################
        #
        # 7. Write TSV file
        #
        ##########################################################
        # lines should be ordered : run, marker, round_down,
        final_optimize_output_df.sort_values(by=['run', 'marker', 'round_down', 'lfn_biosample_replicate: N_ijk/N_jk',
                                                 'biosample', 'replicate'],
                                       ascending=[True, True, True, True, True, True], inplace=True)
        final_optimize_output_df.to_csv(output_file_optimize_lfn, header=True, sep='\t', float_format='%.8f', index=False)
