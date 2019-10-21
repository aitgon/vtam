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
        # Read user known variant information and verify information
        #
        ################################################################################################################

        # variants_optimize_df = pandas.read_csv(input_file_variant_known, sep="\t", header=0, \
        #                                       names=['marker_name', 'run_name', 'biosample_name', 'biosample_type',
        #                                              'variant_id', 'action', 'variant_sequence', 'note'], index_col=False)

        variant_known = VariantKnown(variant_known_tsv, fasta_info_tsv, engine, variant_model, run_model, marker_model, biosample_model, replicate_model)

        # # variants_optimize_df = variant_known.get_variant_known_df()
        #
        # ################################################################################################################
        # #
        # # Control if user variant IDs and sequences are consistent with information in the database
        # #
        # ################################################################################################################
        #
        # variant_control_df = variants_optimize_df[['variant_id', 'variant_sequence']].drop_duplicates()
        # variant_control_df = variant_control_df.loc[~variant_control_df.variant_id.isnull()]
        # with engine.connect() as conn:
        #     for row in variant_control_df.itertuples():
        #         variant_id = row.variant_id
        #         variant_sequence = row.variant_sequence
        #         stmt_select = select([variant_model.__table__.c.id, variant_model.__table__.c.sequence])\
        #             .where(variant_model.__table__.c.id == variant_id)\
        #             .where(variant_model.__table__.c.sequence == variant_sequence)
        #         if conn.execute(stmt_select).first() is None:
        #            Logger.instance().error("Variant {} and its sequence are not coherent with the VTAM database".format(variant_id))
        #            os.mknod(output_file_optimize_pcr_error)
        #            exit()

        ################################################################################################################
        #
        # Read fasta information with current analysis
        #
        ################################################################################################################

        fasta_info = FastaInfo(fasta_info_tsv=fasta_info_tsv, engine=engine)
        fasta_info_records = fasta_info.get_ids_of_run_marker_biosample_replicate(engine, run_model, marker_model, biosample_model, replicate_model)
        fasta_info_df = pandas.DataFrame.from_records(data=fasta_info_records)

        ################################################################################################################
        #
        # Intersect fasta_info and variant_known
        #
        ################################################################################################################


        ################################################################################################################
        #
        # Extract "keep" variants and some columns
        #
        ################################################################################################################
        variant_known_ids_df = variant_known.variant_known_ids_df
        variant_keep_df = variant_known_ids_df.loc[
            variant_known_ids_df['action'] == 'keep', ['variant_id', 'variant_sequence']].drop_duplicates(inplace=False)
        run_marker_biosample_keep_df = variant_known_ids_df.loc[
            variant_known_ids_df['action'] == 'keep', ['run_id', 'marker_id', 'biosample_id']].drop_duplicates(inplace=False)
        variant_keep_df.columns = ['id', 'sequence']

        # import pdb; pdb.set_trace()
        # variants_keep_df = variants_optimize_df[variants_optimize_df["action"].isin(["keep"])]
        # variants_keep_df = variants_keep_df[['marker_name', 'run_name', 'biosample_name', 'variant_id',
        #                                      'variant_sequence']]

        ################################################################################################################
        #
        # For "keep" variants , select run_id, marker_id, variant_id, variant_sequence, biosample_id, replicated_id
        # and read_count (N_ijk)
        #
        ################################################################################################################

        variant_list = []
        with engine.connect() as conn:
            for row in run_marker_biosample_keep_df.itertuples():
                run_id = row.run_id
                marker_id = row.marker_id
                biosample_id = row.biosample_id
                stmt_select = sqlalchemy.select([
                    variant_model.__table__.c.id,
                    variant_model.__table__.c.sequence]) \
                    .where(variant_read_count_model.__table__.c.run_id == run_id) \
                    .where(variant_read_count_model.__table__.c.marker_id == marker_id) \
                    .where(variant_read_count_model.__table__.c.biosample_id == biosample_id) \
                    .where(variant_read_count_model.__table__.c.variant_id == variant_model.__table__.c.id) \
                    .distinct()
                variant_list = variant_list + conn.execute(stmt_select).fetchall()

        variant_df = pandas.DataFrame.from_records(variant_list, columns=['id', 'sequence']).drop_duplicates(inplace=False)

        ################################################################################################################
        #
        # run f10_pcr_error_run_vsearch &  read_count_unexpected_expected_ratio_max
        #
        ################################################################################################################

        vsearch_output_df = f10_pcr_error_run_vsearch(
            variant_db_df=variant_keep_df, variant_usearch_global_unexpected_df=variant_df,
            tmp_dir=this_step_tmp_dir)
        import pdb; pdb.set_trace()

        ################################################################################################################
        #
        #  Create variant information input to run PCRError vsearch
        #
        ################################################################################################################
        variant_read_count_df = variant_read_count_df.rename(columns={'read_count': 'N_ijk'})

        variant_vsearch_db_df = variant_df.loc[variant_df.id.isin(
            variants_keep_df.variant_id.unique().tolist())][['id', 'sequence']].drop_duplicates()

        variant_vsearch_db_df.rename(columns={'variant_id': 'id', 'variant_sequence': 'sequence'}, inplace=True)

        ################################################################################################################
        #
        # run f10_pcr_error_run_vsearch &  read_count_unexpected_expected_ratio_max
        #
        ################################################################################################################

        vsearch_output_df = f10_pcr_error_run_vsearch(
            variant_db_df=variant_vsearch_db_df, variant_usearch_global_unexpected_df=variant_df,
            tmp_dir=this_step_tmp_dir)


        pcr_error_df = f10_get_maximal_pcr_error_value(variant_read_count_df, vsearch_output_df)

        Logger.instance().debug(
            "file: {}; line: {}; pcr error optimize parameter succefully counted : #: {} ".format(__file__,inspect.currentframe().f_lineno, output_file_optimize_pcr_error,'OptimizePCRerror'))
        ##########################################################
        #
        # 7. Write Optimize PCRError to TSV file
        #
        ##########################################################
        pcr_error_df.to_csv(output_file_optimize_pcr_error, header=True, sep='\t', float_format='%.10f', index=False)
