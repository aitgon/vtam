import inspect

import os
import sys

from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
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
            "log_verbosity": "int",
            "log_file": "str",
        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
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
        input_file_variant_known = self.input_file(OptimizePCRerror.__input_file_variant_known)
        input_file_fastainfo = self.input_file(OptimizePCRerror.__input_file_fastainfo)
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
        # Read user known variant information
        #
        ################################################################################################################

        variants_optimize_df = pandas.read_csv(input_file_variant_known, sep="\t", header=0, \
                                              names=['marker_name', 'run_name', 'biosample_name', 'biosample_type',
                                                     'variant_id', 'action', 'variant_sequence', 'note'], index_col=False)

        ################################################################################################################
        #
        # Control if user variant IDs and sequences are consistent with information in the database
        #
        ################################################################################################################

        variant_control_df = variants_optimize_df[['variant_id', 'variant_sequence']].drop_duplicates()
        variant_control_df = variant_control_df.loc[~variant_control_df.variant_id.isnull()]
        with engine.connect() as conn:
            for row in variant_control_df.itertuples():
                variant_id = row.variant_id
                variant_sequence = row.variant_sequence
                stmt_select = select([variant_model.__table__.c.id, variant_model.__table__.c.sequence])\
                    .where(variant_model.__table__.c.id == variant_id)\
                    .where(variant_model.__table__.c.sequence == variant_sequence)
                if conn.execute(stmt_select).first() is None:
                   Logger.instance().error("Variant {} and its sequence are not coherent with the VTAM database".format(variant_id))
                   os.mknod(output_file_optimize_pcr_error)
                   exit()

        ################################################################################################################
        #
        # Extract "keep" variants and some columns
        #
        ################################################################################################################

        variants_keep_df = variants_optimize_df[variants_optimize_df["action"].isin(["keep"])]
        variants_keep_df = variants_keep_df[['marker_name', 'run_name', 'biosample_name', 'variant_id',
                                             'variant_sequence']]

        ################################################################################################################
        #
        # For "keep" variants , select run_id, marker_id, variant_id, variant_sequence, biosample_id, replicated_id
        # and read_count (N_ijk)
        #
        ################################################################################################################

        variant_list = []
        with engine.connect() as conn:
            for row in variants_keep_df.itertuples():
                run_name = row.run_name
                marker_name = row.marker_name
                biosample_name = row.biosample_name
                stmt_select = select([
                    variant_model.__table__.c.id,
                    variant_model.__table__.c.sequence, ]) \
                    .where(run_model.__table__.c.name == run_name) \
                    .where(variant_read_count_model.__table__.c.run_id == run_model.__table__.c.id) \
                    .where(marker_model.__table__.c.name == marker_name) \
                    .where(variant_read_count_model.__table__.c.marker_id == marker_model.__table__.c.id) \
                    .where(biosample_model.__table__.c.name == biosample_name) \
                    .where(variant_read_count_model.__table__.c.biosample_id == biosample_model.__table__.c.id) \
                    .where(variant_read_count_model.__table__.c.variant_id == variant_model.__table__.c.id) \
                    .distinct()
                variant_list = variant_list + conn.execute(stmt_select).fetchall()

        variant_df = pandas.DataFrame.from_records(variant_list, columns=['id', 'sequence'])

        ##########################################################
        #
        # Read fastainfo to get run_id, marker_id, biosample_id, replicate_id for current analysis
        #
        ##########################################################
        fastainfo_df = pandas.read_csv(input_file_fastainfo, sep="\t", header=0, \
                                       names=['tag_forward', 'primer_forward', 'tag_reverse', 'primer_reverse',
                                              'marker_name', 'biosample_name', \
                                              'replicate_name', 'run_name', 'fastq_fwd', 'fastq_rev', 'fasta'])
        sample_instance_list = []
        for row in fastainfo_df.itertuples():
            marker_name = row.marker_name
            run_name = row.run_name
            biosample_name = row.biosample_name
            replicate_name = row.replicate_name
            with engine.connect() as conn:
                # get run_id ###########
                stmt_select_run_id = select([run_model.__table__.c.id]).where(run_model.__table__.c.name == run_name)
                run_id = conn.execute(stmt_select_run_id).first()[0]
                # get marker_id ###########
                stmt_select_marker_id = select([marker_model.__table__.c.id]).where(
                    marker_model.__table__.c.name == marker_name)
                marker_id = conn.execute(stmt_select_marker_id).first()[0]
                # get biosample_id ###########
                stmt_select_biosample_id = select([biosample_model.__table__.c.id]).where(
                    biosample_model.__table__.c.name == biosample_name)
                biosample_id = conn.execute(stmt_select_biosample_id).first()[0]
                # get replicate_id ###########
                stmt_select_replicate_id = select([replicate_model.__table__.c.id]).where(
                    replicate_model.__table__.c.name == replicate_name)
                replicate_id = conn.execute(stmt_select_replicate_id).first()[0]
                # add this sample_instance ###########
                sample_instance_list.append({'run_id': run_id, 'marker_id': marker_id, 'biosample_id': biosample_id,
                                             'replicate_id': replicate_id})

        ##########################################################
        #
        # Select marker/run/biosample/replicate from variant_read_count_model
        #
        ##########################################################

        variant_read_count_model_table = variant_read_count_model.__table__

        variant_read_count_list = []
        for sample_instance in sample_instance_list:
            run_id = sample_instance['run_id']
            marker_id = sample_instance['marker_id']
            biosample_id = sample_instance['biosample_id']
            replicate_id = sample_instance['replicate_id']
            stmt_select = select([variant_read_count_model_table.c.run_id,
                                  variant_read_count_model_table.c.marker_id,
                                  variant_read_count_model_table.c.biosample_id,
                                  variant_read_count_model_table.c.replicate_id,
                                  variant_read_count_model_table.c.variant_id,
                                  variant_read_count_model_table.c.read_count]).distinct() \
                .where(variant_read_count_model.__table__.c.run_id == run_id) \
                .where(variant_read_count_model.__table__.c.marker_id == marker_id) \
                .where(variant_read_count_model.__table__.c.biosample_id == biosample_id) \
                .where(variant_read_count_model.__table__.c.replicate_id == replicate_id)
            with engine.connect() as conn:
                for row2 in conn.execute(stmt_select).fetchall():
                    variant_read_count_list.append(row2)
        #
        variant_read_count_df = pandas.DataFrame.from_records(variant_read_count_list,
                                                              columns=['run_id', 'marker_id', 'variant_id',
                                                                       'biosample_id', 'replicate_id', 'read_count'])

        # Exit if no variants for analysis
        try:
            assert variant_read_count_df.shape[0] > 0
        except AssertionError:
            sys.stderr.write("Error: No variants available for this filter: {}".format(os.path.basename(__file__)))
            sys.exit(1)

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
            variant_db_df=variant_vsearch_db_df, variant_usearch_global_df=variant_df,
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
