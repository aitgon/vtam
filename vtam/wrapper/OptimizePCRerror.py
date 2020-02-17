import os
import pathlib

import sqlalchemy
from wopmars.models.ToolWrapper import ToolWrapper

from vtam.utils.FilterPCRerrorRunner import FilterPCRerrorRunner
from vtam.utils.SampleInformationUtils import SampleInformationUtils, FastaInformationTSV
from vtam.utils.VariantKnown import VariantKnown
import pandas

from vtam.utils.PathManager import PathManager


class OptimizePCRerror(ToolWrapper):

    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.OptimizePCRerror"
    }

    # Input file
    __input_file_fastainfo = "fastainfo"
    __input_file_variant_known = "variant_known"
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
            OptimizePCRerror.__input_file_fastainfo,
            OptimizePCRerror.__input_file_variant_known,
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

        this_temp_dir = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__))
        pathlib.Path(this_temp_dir).mkdir(exist_ok=True)

        ################################################################################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ################################################################################################################

        # Input file paths
        variant_known_tsv = self.input_file(OptimizePCRerror.__input_file_variant_known)
        fasta_info_tsv = self.input_file(OptimizePCRerror.__input_file_fastainfo)
        #
        # Input models
        run_model = self.input_table(OptimizePCRerror.__input_table_run)
        marker_model = self.input_table(OptimizePCRerror.__input_table_marker)
        variant_model = self.input_table(OptimizePCRerror.__input_table_variant)
        biosample_model = self.input_table(OptimizePCRerror.__input_table_biosample)
        variant_read_count_model = self.input_table(OptimizePCRerror.__input_table_variant_read_count)
        #
        # Output file paths
        output_file_optimize_pcr_error = self.output_file(OptimizePCRerror.__output_file_optimize_pcr_error)

        ################################################################################################################
        #
        # Read user known variant information and verify consistency with DB
        #
        ################################################################################################################

        variant_known = VariantKnown(variant_known_tsv, fasta_info_tsv, engine)
        variant_known_ids_df = variant_known.variant_known_ids_df

        ################################################################################################################
        #
        # Get run_marker_biosample IDs marked as mock
        #
        ################################################################################################################

        run_marker_biosample_mock_df = variant_known_ids_df.loc[
            variant_known_ids_df.biosample_type == 'mock', ['run_id', 'marker_id', 'biosample_id']]
        run_marker_biosample_mock_df.drop_duplicates(inplace=True)

        ##########################################################
        #
        # variant_read_count_df
        #
        ##########################################################

        fasta_info_tsv = FastaInformationTSV(fasta_info_tsv=fasta_info_tsv, engine=engine)
        sample_information_df_analyzer = SampleInformationUtils(engine, fasta_info_tsv.sample_information_df)
        variant_read_count_df = sample_information_df_analyzer.get_variant_read_count_df(variant_read_count_like_model=variant_read_count_model)

        variant_df = sample_information_df_analyzer.get_variant_df(variant_read_count_like_model=variant_read_count_model, variant_model=variant_model)

        ################################################################################################################
        #
        # Get variant_read_count_expected, variant_read_count_unexpected and variant_df
        #
        ################################################################################################################

        run_marker_biosample_variant_keep_df = variant_known.get_keep_run_marker_biosample_variant_df()

        variant_delete_df, variant_delete_mock_df, variant_delete_negative_df, variant_delete_real_df\
            = variant_known.get_delete_run_marker_biosample_variant_df(variant_read_count_df=variant_read_count_df)

        ################################################################################################################
        #
        # Run per biosample_id
        #
        ################################################################################################################

        pcr_error_final_df = pandas.DataFrame()

        for biosample_id in run_marker_biosample_mock_df.biosample_id.unique().tolist():

            run_marker_biosample_variant_keep_per_biosample_df = run_marker_biosample_variant_keep_df.loc[
                run_marker_biosample_variant_keep_df.biosample_id==biosample_id]
            variant_expected_per_biosample_df = variant_df.loc[
                variant_df.index.isin(run_marker_biosample_variant_keep_per_biosample_df.variant_id.unique().tolist())].drop_duplicates()

            variant_delete_mock_per_biosample_df = variant_delete_mock_df.loc[variant_delete_mock_df.biosample_id == biosample_id]
            variant_unexpected_per_biosample_df = variant_df.loc[
                variant_df.index.isin(variant_delete_mock_per_biosample_df.variant_id.unique().tolist())].drop_duplicates()

            variant_read_count_per_biosample_df = variant_read_count_df.loc[variant_read_count_df.biosample_id == biosample_id]

            this_step_tmp_per_biosample_dir = os.path.join(this_temp_dir, str(biosample_id))
            pathlib.Path(this_step_tmp_per_biosample_dir).mkdir(exist_ok=True)

            ##########################################################
            #
            # Run vsearch and get alignement df
            #
            ##########################################################

            filter_pcr_error_runner = FilterPCRerrorRunner(variant_expected_df=variant_expected_per_biosample_df,
                                                           variant_unexpected_df=variant_unexpected_per_biosample_df,
                                                           variant_read_count_df=variant_read_count_per_biosample_df,
                                                           tmp_dir=this_step_tmp_per_biosample_dir)

            pcr_error_df = filter_pcr_error_runner.get_variant_unexpected_to_expected_ratio_df()

            pcr_error_final_df = pandas.concat([pcr_error_final_df, pcr_error_df], axis=0)

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
        pcr_error_final_df = pcr_error_final_df.merge(run_id_to_name_df, on='run_id')

        marker_id_to_name_df = pandas.DataFrame.from_records(data=marker_id_to_name, columns=['marker_id', 'marker_name'])
        pcr_error_final_df = pcr_error_final_df.merge(marker_id_to_name_df, on='marker_id')

        biosample_id_to_name_df = pandas.DataFrame.from_records(data=biosample_id_to_name, columns=['biosample_id', 'biosample_name'])
        pcr_error_final_df = pcr_error_final_df.merge(biosample_id_to_name_df, on='biosample_id')

        pcr_error_final_df = pcr_error_final_df[['run_name', 'marker_name', 'biosample_name', 'variant_id_expected',
                'N_ij_expected', 'variant_id_unexpected', 'N_ij_unexpected', 'N_ij_unexpected_to_expected_ratio']]

        ##########################################################
        #
        # Write Optimize PCRError to TSV file
        #
        ##########################################################

        pcr_error_final_df.sort_values(by=['N_ij_unexpected_to_expected_ratio', 'variant_id_expected', 'variant_id_unexpected'],
                                       ascending=[False, True, True], inplace=True)
        pcr_error_final_df.to_csv(output_file_optimize_pcr_error, header=True, sep='\t', float_format='%.10f', index=False)

