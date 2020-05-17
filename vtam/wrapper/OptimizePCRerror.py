import os
import pathlib

import sqlalchemy
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

        #######################################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        #######################################################################

        # Input file paths
        known_occurrences_tsv = self.input_file(
            OptimizePCRerror.__input_file_known_occurrences)
        fasta_info_tsv_path = self.input_file(
            OptimizePCRerror.__input_file_readinfo)
        #
        # Output file paths
        output_file_optimize_pcr_error = self.output_file(
            OptimizePCRerror.__output_file_optimize_pcr_error)

        #######################################################################
        #
        # Group and run this genetic_code by run/marker combination
        # Loop by run/marker
        #
        #######################################################################

        final_pcr_error_df = pandas.DataFrame()

        known_occurrences_df = pandas.read_csv(
            known_occurrences_tsv, sep="\t", header=0)
        known_occurrences_df.columns = known_occurrences_df.columns.str.lower()
        vknown_grouped = known_occurrences_df.groupby(by=['run', 'marker'])
        for vknown_grouped_key in vknown_grouped.groups:
            known_occurrences_i_df = known_occurrences_df.loc[
                vknown_grouped.groups[vknown_grouped_key], :]

            ###################################################################
            #
            # Read user known variant information and verify consistency with DB
            #
            ###################################################################

            known_occurrences = KnownOccurrences(
                known_occurrences_i_df, fasta_info_tsv_path, engine)
            known_occurrences_ids_df = known_occurrences.known_occurrences_ids_df

            ###################################################################
            #
            # Get run_marker_biosample IDs marked as mock
            #
            ###################################################################

            run_marker_biosample_mock_df = known_occurrences_ids_df.loc[
                known_occurrences_ids_df.biosample_type == 'mock', ['run_id', 'marker_id', 'biosample_id']]
            run_marker_biosample_mock_df.drop_duplicates(inplace=True)

            ##########################################################
            #
            # variant_read_count_input_df
            #
            ##########################################################

            # fasta_info_tsv_obj = FastaInformationTSV(fasta_info_tsv=fasta_info_tsv_path, engine=engine)
            # sample_information_df_analyzer = SampleInformationUtils(engine, fasta_info_tsv_obj.sample_information_df)
            # variant_read_count_df = sample_information_df_analyzer.get_variant_read_count_df(variant_read_count_like_model=VariantReadCount)
            sample_info_tsv_obj = SampleInformationFile(
                tsv_path=fasta_info_tsv_path)
            variant_read_count_df = sample_info_tsv_obj.get_variant_read_count_df(
                VariantReadCount, engine=engine)

            # variant_df = sample_information_df_analyzer.get_variant_df(variant_read_count_like_model=VariantReadCount, variant_model=Variant)
            variant_df = pandas.read_sql(
                sqlalchemy.select(
                    [
                        Variant.__table__.c.id,
                        Variant.__table__.c.sequence]),
                con=engine.connect(),
                index_col='id')

            ###################################################################
            #
            # Get variant_read_count_expected, variant_read_count_unexpected and variant_df
            #
            ###################################################################

            run_marker_biosample_variant_keep_df = known_occurrences.get_keep_run_marker_biosample_variant_df()

            variant_delete_df, variant_delete_mock_df, variant_delete_negative_df, variant_delete_real_df\
                = known_occurrences.get_delete_run_marker_biosample_variant_df(variant_read_count_df=variant_read_count_df)

            ###################################################################
            #
            # Run per biosample_id
            #
            ###################################################################

            pcr_error_per_run_marker_df = pandas.DataFrame()

            for biosample_id in run_marker_biosample_mock_df.biosample_id.unique().tolist():

                run_marker_biosample_variant_keep_per_biosample_df = run_marker_biosample_variant_keep_df.loc[
                    run_marker_biosample_variant_keep_df.biosample_id == biosample_id]
                variant_expected_per_biosample_df = variant_df.loc[variant_df.index.isin(
                    run_marker_biosample_variant_keep_per_biosample_df.variant_id.unique().tolist())].drop_duplicates()

                variant_delete_mock_per_biosample_df = variant_delete_mock_df.loc[
                    variant_delete_mock_df.biosample_id == biosample_id]
                variant_unexpected_per_biosample_df = variant_df.loc[variant_df.index.isin(
                    variant_delete_mock_per_biosample_df.variant_id.unique().tolist())].drop_duplicates()

                variant_read_count_per_biosample_df = variant_read_count_df.loc[
                    variant_read_count_df.biosample_id == biosample_id]

                this_step_tmp_per_biosample_dir = os.path.join(
                    this_temp_dir, str(biosample_id))
                pathlib.Path(this_step_tmp_per_biosample_dir).mkdir(
                    exist_ok=True)

                ###############################################################
                #
                # Run vsearch and get alignement variant_read_count_input_df
                #
                ###############################################################

                filter_pcr_error_runner = FilterPCRerrorRunner(
                    variant_expected_df=variant_expected_per_biosample_df,
                    variant_unexpected_df=variant_unexpected_per_biosample_df,
                    variant_read_count_df=variant_read_count_per_biosample_df,
                    tmp_dir=this_step_tmp_per_biosample_dir)

                pcr_error_df = filter_pcr_error_runner.get_variant_unexpected_to_expected_ratio_df()

                pcr_error_per_run_marker_df = pandas.concat(
                    [pcr_error_per_run_marker_df, pcr_error_df], axis=0)

            ###################################################################
            #
            # Convert run_id, marker_id and biosample_id to their names
            #
            ###################################################################

            with engine.connect() as conn:
                run_id_to_name = conn.execute(sqlalchemy.select(
                    [Run.__table__.c.id, Run.__table__.c.name])).fetchall()
                marker_id_to_name = conn.execute(sqlalchemy.select(
                    [Marker.__table__.c.id, Marker.__table__.c.name])).fetchall()
                biosample_id_to_name = conn.execute(sqlalchemy.select(
                    [Biosample.__table__.c.id, Biosample.__table__.c.name])).fetchall()

            run_id_to_name_df = pandas.DataFrame.from_records(
                data=run_id_to_name, columns=['run_id', 'run_name'])
            pcr_error_per_run_marker_df = pcr_error_per_run_marker_df.merge(
                run_id_to_name_df, on='run_id')

            marker_id_to_name_df = pandas.DataFrame.from_records(
                data=marker_id_to_name, columns=['marker_id', 'marker_name'])
            pcr_error_per_run_marker_df = pcr_error_per_run_marker_df.merge(
                marker_id_to_name_df, on='marker_id')

            biosample_id_to_name_df = pandas.DataFrame.from_records(
                data=biosample_id_to_name, columns=['biosample_id', 'biosample_name'])
            pcr_error_per_run_marker_df = pcr_error_per_run_marker_df.merge(
                biosample_id_to_name_df, on='biosample_id')

            pcr_error_per_run_marker_df = pcr_error_per_run_marker_df[['run_name',
                                                                       'marker_name',
                                                                       'biosample_name',
                                                                       'variant_id_expected',
                                                                       'N_ij_expected',
                                                                       'variant_id_unexpected',
                                                                       'N_ij_unexpected',
                                                                       'N_ij_unexpected_to_expected_ratio']]

            final_pcr_error_df = pandas.concat(
                [final_pcr_error_df, pcr_error_per_run_marker_df], axis=0)

        #######################################################################
        #
        # Add sequence
        #
        #######################################################################

        final_pcr_error_df.rename({'marker_name': 'marker',
                                   'run_name': 'run',
                                   'biosample_name': 'biosample'},
                                  axis=1,
                                  inplace=True)
        final_pcr_error_df['variant_seq_unexpected'] = ''
        with engine.connect() as conn:
            for variant_id_unexpected in final_pcr_error_df.variant_id_unexpected.unique():
                variant_id_unexpected = int(variant_id_unexpected)
                variant_sequence_row = conn.execute(sqlalchemy.select([Variant.__table__.c.sequence]).where(
                    Variant.__table__.c.id == variant_id_unexpected)).first()
                if not (variant_sequence_row is None):
                    variant_sequence = variant_sequence_row[0]
                    final_pcr_error_df.loc[(final_pcr_error_df.variant_id_unexpected ==
                                            variant_id_unexpected).values, 'variant_seq_unexpected'] = variant_sequence

        final_pcr_error_df.rename({'marker_name': 'marker',
                                   'run_name': 'run',
                                   'biosample_name': 'biosample'},
                                  axis=1,
                                  inplace=True)
        final_pcr_error_df['variant_seq_expected'] = ''
        with engine.connect() as conn:
            for variant_id_expected in final_pcr_error_df.variant_id_expected.unique():
                variant_id_expected = int(variant_id_expected)
                variant_sequence_row = conn.execute(sqlalchemy.select([Variant.__table__.c.sequence]).where(
                    Variant.__table__.c.id == variant_id_expected)).first()
                if not (variant_sequence_row is None):
                    variant_sequence = variant_sequence_row[0]
                    final_pcr_error_df.loc[(final_pcr_error_df.variant_id_expected ==
                                            variant_id_expected).values, 'variant_seq_expected'] = variant_sequence

        #######################################################################
        #
        # Write Optimize PCRError to TSV file
        #
        #######################################################################

        # lines should be ordered: by run, marker,
        final_pcr_error_df.sort_values(by=['run',
                                           'marker',
                                           'N_ij_unexpected_to_expected_ratio',
                                           'variant_id_expected',
                                           'variant_id_unexpected'],
                                       ascending=[True,
                                                  True,
                                                  False,
                                                  True,
                                                  True],
                                       inplace=True)
        final_pcr_error_df.to_csv(
            output_file_optimize_pcr_error,
            header=True,
            sep='\t',
            float_format='%.10f',
            index=False)
