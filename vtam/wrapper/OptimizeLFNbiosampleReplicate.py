import os
import sys

from wopmars.framework.database.tables.ToolWrapper import ToolWrapper

from vtam.utils.OptionManager import OptionManager
from sqlalchemy import select
import pandas

from vtam.utils.Logger import Logger


class OptimizeLFNbiosampleReplicate(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.OptimizeLFNbiosampleReplicate"
    }

    # Input file
    __input_file_fastainfo = "fastainfo"
    __input_file_variant_known = "variant_known"
    # Input table
    __input_table_run = "Run"
    __input_table_marker = "Marker"
    __input_table_biosample = "Biosample"
    __input_table_replicate = "Replicate"
    __input_table_variant = "Variant"
    __input_table_variant_read_count = "VariantReadCount"
    # Output file
    __output_file_optimize_lfn_biosample_replicate = "optimize_lfn_biosample_replicate"


    def specify_input_file(self):
        return[
            OptimizeLFNbiosampleReplicate.__input_file_fastainfo,
            OptimizeLFNbiosampleReplicate.__input_file_variant_known,
        ]

    def specify_input_table(self):
        return [
            OptimizeLFNbiosampleReplicate.__input_table_marker,
            OptimizeLFNbiosampleReplicate.__input_table_run,
            OptimizeLFNbiosampleReplicate.__input_table_biosample,
            OptimizeLFNbiosampleReplicate.__input_table_replicate,
            OptimizeLFNbiosampleReplicate.__input_table_variant,
            OptimizeLFNbiosampleReplicate.__input_table_variant_read_count,
        ]


    def specify_output_file(self):
        return [
            OptimizeLFNbiosampleReplicate.__output_file_optimize_lfn_biosample_replicate,
        ]

    def specify_params(self):
        return {
            "log_verbosity": "int",
            "log_file": "str"
        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        OptionManager.instance()['log_verbosity'] = int(self.option("log_verbosity"))
        OptionManager.instance()['log_file'] = str(self.option("log_file"))

        ##########################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ##########################################################
        #
        # Input file path
        input_file_variant_known = self.input_file(OptimizeLFNbiosampleReplicate.__input_file_variant_known)
        input_file_fastainfo = self.input_file(OptimizeLFNbiosampleReplicate.__input_file_fastainfo)
        #
        # Input table models
        run_model = self.input_table(OptimizeLFNbiosampleReplicate.__input_table_run)
        marker_model = self.input_table(OptimizeLFNbiosampleReplicate.__input_table_marker)
        biosample_model = self.input_table(OptimizeLFNbiosampleReplicate.__input_table_biosample)
        replicate_model = self.input_table(OptimizeLFNbiosampleReplicate.__input_table_replicate)
        variant_model = self.input_table(OptimizeLFNbiosampleReplicate.__input_table_variant)
        variant_read_count_model = self.input_table(OptimizeLFNbiosampleReplicate.__input_table_variant_read_count)
        #
        # Output file path
        output_file_optimize_lfn = self.output_file(OptimizeLFNbiosampleReplicate.__output_file_optimize_lfn_biosample_replicate)

        ##########################################################
        #
        # 1. Read variants_optimize to get run_id, marker_id, biosample_id, variant_id for current analysis
        #
        ##########################################################
        # positive_variant_df = pandas.read_csv(input_file_variant_known, sep="\t", header=0,\
        #     names=['marker_name', 'run_name', 'biosample_name', 'variant_id', 'variant_sequence'], index_col=False)
        variants_optimize_df = pandas.read_csv(input_file_variant_known, sep="\t", header=0, \
                                              names=['marker_name', 'run_name', 'biosample_name', 'biosample_type',
                                                     'variant_id', 'action', 'variant_sequence', 'note'], index_col=False)

        ########################
        #
        # Control if user variants and sequence are consistent in the database
        #
        ########################
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
                   os.mknod(output_file_optimize_lfn)
                   exit()


        ##########################################################
        #
        # Extract some columns and "keep" variants
        #
        ##########################################################
        variants_keep_tolerate_df = variants_optimize_df[variants_optimize_df["action"].isin(["keep", "tolerate"])]
        variants_keep_tolerate_df = variants_keep_tolerate_df[['marker_name', 'run_name', 'biosample_name', 'variant_id', 'variant_sequence']]

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
                                                              columns=['run_id', 'marker_id', 'biosample_id',
                                                                       'replicate_id', 'variant_id', 'read_count'])

        # Exit if no variants for analysis
        try:
            assert variant_read_count_df.shape[0] > 0
        except AssertionError:
            sys.stderr.write("Error: No variants available for this filter: {}".format(os.path.basename(__file__)))
            sys.exit(1)


        ##########################################################
        #
        # 4. Compute ratio per_sum_biosample_replicate: N_ijk / N_jk
        #
        ##########################################################
        variant_read_count_df = variant_read_count_df.rename(columns={'read_count': 'N_ijk'})
        aggregate_df = variant_read_count_df[['run_id', 'marker_id', 'biosample_id', 'replicate_id', 'N_ijk']].groupby(
            by=['run_id', 'marker_id', 'biosample_id', 'replicate_id']).sum().reset_index()
        aggregate_df = aggregate_df.rename(columns={'N_ijk': 'N_jk'})
        variant_read_count_df = variant_read_count_df.merge(aggregate_df, on=['run_id', 'marker_id', 'biosample_id', 'replicate_id'])
        variant_read_count_df['lfn_biosample_replicate: N_ijk/N_jk'] = variant_read_count_df['N_ijk'] / variant_read_count_df['N_jk']

        ##########################################################
        #
        # 5.Sort and extract the lowest value of the ration with 4 digit
        #
        ##########################################################
        variant_read_count_df=variant_read_count_df.sort_values('lfn_biosample_replicate: N_ijk/N_jk', ascending=True)
        # Make round.inf with 4 decimals
        round_inf_4_decimals = lambda x: int(x * 10 ** 4) / 10 ** 4
        variant_read_count_df['round_inf'] = variant_read_count_df['lfn_biosample_replicate: N_ijk/N_jk'].apply(round_inf_4_decimals)

        ##########################################################
        #
        # 7. Write TSV file
        #
        ##########################################################
        variant_read_count_df.to_csv(output_file_optimize_lfn, header=True, sep='\t', float_format='%.10f', index=False)

