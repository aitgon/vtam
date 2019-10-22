import inspect

import os
import sys

import sqlalchemy
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper

from vtam.utils.FastaInfo import FastaInfo
from vtam.utils.VariantKnown import VariantKnown
from vtam.wrapper.FilterPCRError import f10_get_maximal_pcr_error_value, f10_pcr_error_run_vsearch

from sqlalchemy import select
import pandas

from vtam.utils.PathManager import PathManager

from vtam.utils.Logger import Logger
from vtam.utils.OptionManager import OptionManager



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
    __input_table_replicate = "Replicate"
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
            OptimizePCRerror.__input_table_replicate,
            OptimizePCRerror.__input_table_variant_read_count,
        ]


    def specify_output_file(self):
        return [
            OptimizePCRerror.__output_file_optimize_pcr_error,
        ]

    def specify_params(self):
        return {
            "foo": "int",
            "log_verbosity": "int",
            "log_file": "str",
        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        if not self.option("log_verbosity") is None:
            OptionManager.instance()['log_verbosity'] = int(self.option("log_verbosity"))
            OptionManager.instance()['log_file'] = str(self.option("log_file"))
        this_step_tmp_dir = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__))
        PathManager.mkdir_p(this_step_tmp_dir)

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
        replicate_model = self.input_table(OptimizePCRerror.__input_table_replicate)
        variant_read_count_model = self.input_table(OptimizePCRerror.__input_table_variant_read_count)
        #
        # Output file paths
        output_file_optimize_pcr_error = self.output_file(OptimizePCRerror.__output_file_optimize_pcr_error)

        ################################################################################################################
        #
        # Read fasta information with current analysis
        #
        ################################################################################################################

        fasta_info = FastaInfo(fasta_info_tsv=fasta_info_tsv, engine=engine)

        ################################################################################################################
        #
        # Read user known variant information and verify information
        #
        ################################################################################################################

        variant_known = VariantKnown(variant_known_tsv, fasta_info_tsv, engine, variant_model, run_model, marker_model,
                                     biosample_model, replicate_model)
        variant_known_ids_df = variant_known.variant_known_ids_df

        # aggregate per mock biosample
        run_marker_biosample_aggregated_by_mock_df = variant_known_ids_df.loc[
            variant_known_ids_df.biosample_type == 'mock', ['run_id', 'marker_id', 'biosample_id']]
        run_marker_biosample_aggregated_by_mock_df.drop_duplicates(inplace=True)

        ################################################################################################################
        #
        # Run per biosample mock
        #
        ################################################################################################################
        pcr_error_final_df = pandas.DataFrame()
        for per_biosample_mock in run_marker_biosample_aggregated_by_mock_df.itertuples():
            per_biosample_mock_run_id = per_biosample_mock.run_id
            per_biosample_mock_marker_id = per_biosample_mock.marker_id
            per_biosample_mock_biosample_id = per_biosample_mock.biosample_id

            ################################################################################################################
            #
            # Extract "keep" and "tolerate" variants and some columns
            #
            ################################################################################################################

            variant_keep_tolerate_df = variant_known_ids_df.loc[
                (variant_known_ids_df['action'].isin(['keep', 'tolerate'])) &
                (variant_known_ids_df['run_id'] == per_biosample_mock_run_id) &
                (variant_known_ids_df['marker_id'] == per_biosample_mock_marker_id) &
                (variant_known_ids_df['biosample_id'] == per_biosample_mock_biosample_id),
                ['variant_id', 'variant_sequence']]
            variant_keep_tolerate_df.drop_duplicates(inplace=True)
            variant_keep_tolerate_df.variant_id = variant_keep_tolerate_df.variant_id.astype('int')

            variant_keep_tolerate_df.columns = ['id', 'sequence']

            ################################################################################################################
            #
            # For run, marker and mock biosamples, create the variant df with id and sequence
            #
            ################################################################################################################

            variant_list = []
            with engine.connect() as conn:
                stmt_select = sqlalchemy.select([
                    variant_model.__table__.c.id,
                    variant_model.__table__.c.sequence]) \
                    .where(variant_read_count_model.__table__.c.run_id == per_biosample_mock_run_id) \
                    .where(variant_read_count_model.__table__.c.marker_id == per_biosample_mock_marker_id) \
                    .where(variant_read_count_model.__table__.c.biosample_id == per_biosample_mock_biosample_id) \
                    .where(variant_read_count_model.__table__.c.variant_id == variant_model.__table__.c.id) \
                    .distinct()
                variant_list = variant_list + conn.execute(stmt_select).fetchall()

            variant_df = pandas.DataFrame.from_records(variant_list, columns=['id', 'sequence']).drop_duplicates(inplace=False)

            ################################################################################################################
            #
            # For run, marker and mock biosamples, create the variant_read_count df
            #
            ################################################################################################################

            variant_read_count_list = []
            with engine.connect() as conn:
                stmt_select = sqlalchemy.select([
                    variant_read_count_model.__table__.c.run_id,
                    variant_read_count_model.__table__.c.marker_id,
                    variant_read_count_model.__table__.c.biosample_id,
                    variant_read_count_model.__table__.c.replicate_id,
                    variant_read_count_model.__table__.c.variant_id,
                    variant_read_count_model.__table__.c.read_count,
                ]) \
                    .where(variant_read_count_model.__table__.c.run_id == per_biosample_mock_run_id) \
                    .where(variant_read_count_model.__table__.c.marker_id == per_biosample_mock_marker_id) \
                    .where(variant_read_count_model.__table__.c.biosample_id == per_biosample_mock_biosample_id)\
                    .distinct()
                variant_read_count_list = variant_read_count_list + conn.execute(stmt_select).fetchall()

            variant_read_count_df = pandas.DataFrame.from_records(variant_read_count_list, columns=['run_id', 'marker_id',
                                                            'biosample_id', 'replicate_id', 'variant_id', 'read_count'])
            variant_read_count_df.drop_duplicates(inplace=True)

            ################################################################################################################
            #
            # run f10_pcr_error_run_vsearch &  read_count_unexpected_expected_ratio_max
            #
            ################################################################################################################

            vsearch_output_df = f10_pcr_error_run_vsearch(
                variant_db_df=variant_keep_tolerate_df, variant_usearch_global_unexpected_df=variant_df,
                tmp_dir=this_step_tmp_dir)

            variant_read_count_df = variant_read_count_df.rename(columns={'read_count': 'N_ijk'})
            pcr_error_df = f10_get_maximal_pcr_error_value(variant_read_count_df, vsearch_output_df)

            pcr_error_final_df = pandas.concat([pcr_error_final_df, pcr_error_df], axis=0)

        ##########################################################
        #
        # 7. Write Optimize PCRError to TSV file
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
                'N_ijk_expected', 'variant_id_unexpected', 'N_ijk_unexpected', 'N_ijk_unexpected_expected_ratio']]

        pcr_error_final_df.sort_values(by=['N_ijk_unexpected_expected_ratio', 'variant_id_expected', 'variant_id_unexpected'],
                                       ascending=[False, True, True], inplace=True)
        pcr_error_final_df.to_csv(output_file_optimize_pcr_error, header=True, sep='\t', float_format='%.10f', index=False)

