import os
import sys

import sqlalchemy
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper

from vtam.utils.OptionManager import OptionManager
from sqlalchemy import select
import pandas

from vtam.utils.Logger import Logger
from vtam.utils.VariantKnown import VariantKnown


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
            "foo": "int",
            "log_verbosity": "int",
            "log_file": "str"
        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        if not self.option("log_verbosity") is None:
            OptionManager.instance()['log_verbosity'] = int(self.option("log_verbosity"))
            OptionManager.instance()['log_file'] = str(self.option("log_file"))

        ##########################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ##########################################################

        # Input file output
        variant_known_tsv = self.input_file(OptimizeLFNbiosampleReplicate.__input_file_variant_known)
        fasta_info_tsv = self.input_file(OptimizeLFNbiosampleReplicate.__input_file_fastainfo)
        #
        # Input table models
        run_model = self.input_table(OptimizeLFNbiosampleReplicate.__input_table_run)
        marker_model = self.input_table(OptimizeLFNbiosampleReplicate.__input_table_marker)
        biosample_model = self.input_table(OptimizeLFNbiosampleReplicate.__input_table_biosample)
        replicate_model = self.input_table(OptimizeLFNbiosampleReplicate.__input_table_replicate)
        variant_model = self.input_table(OptimizeLFNbiosampleReplicate.__input_table_variant)
        variant_read_count_model = self.input_table(OptimizeLFNbiosampleReplicate.__input_table_variant_read_count)
        #
        # Output file output
        output_file_optimize_lfn = self.output_file(OptimizeLFNbiosampleReplicate.__output_file_optimize_lfn_biosample_replicate)

        ################################################################################################################
        #
        # Read user known variant information and verify information
        #
        ################################################################################################################

        variant_known = VariantKnown(variant_known_tsv, fasta_info_tsv, engine, variant_model, run_model, marker_model,
                                     biosample_model, replicate_model)
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
        # For run, marker and mock biosamples, create the variant_read_count df
        #
        ################################################################################################################

        variant_read_count_list = []
        with engine.connect() as conn:
            for row in run_marker_biosample_mock_df.itertuples():
                run_id = row.run_id
                marker_id = row.marker_id
                biosample_id = row.biosample_id
                stmt_select = sqlalchemy.select([
                    variant_read_count_model.__table__.c.run_id,
                    variant_read_count_model.__table__.c.marker_id,
                    variant_read_count_model.__table__.c.biosample_id,
                    variant_read_count_model.__table__.c.replicate_id,
                    variant_read_count_model.__table__.c.variant_id,
                    variant_read_count_model.__table__.c.read_count,
                ]) \
                    .where(variant_read_count_model.__table__.c.run_id == run_id) \
                    .where(variant_read_count_model.__table__.c.marker_id == marker_id) \
                    .where(variant_read_count_model.__table__.c.biosample_id == biosample_id)\
                    .distinct()
                variant_read_count_list = variant_read_count_list + conn.execute(stmt_select).fetchall()

        variant_read_count_df = pandas.DataFrame.from_records(variant_read_count_list, columns=['run_id', 'marker_id',
                                                        'biosample_id', 'replicate_id', 'variant_id', 'read_count'])
        variant_read_count_df.drop_duplicates(inplace=True)


        ##########################################################
        #
        # 4. Compute ratio per_sum_biosample_replicate: N_ijk / N_jk
        #
        ##########################################################

        output_df = variant_read_count_df.rename(columns={'read_count': 'N_ijk'})
        aggregate_df = output_df[['run_id', 'marker_id', 'biosample_id', 'replicate_id', 'N_ijk']].groupby(
            by=['run_id', 'marker_id', 'biosample_id', 'replicate_id']).sum().reset_index()
        aggregate_df = aggregate_df.rename(columns={'N_ijk': 'N_jk'})
        output_df = output_df.merge(aggregate_df, on=['run_id', 'marker_id', 'biosample_id', 'replicate_id'])
        output_df['lfn_biosample_replicate: N_ijk/N_jk'] = output_df['N_ijk'] / output_df['N_jk']

        ##########################################################
        #
        # 5.Sort and extract the lowest value of the ration with 4 digit
        #
        ##########################################################

        output_df=output_df.sort_values('lfn_biosample_replicate: N_ijk/N_jk', ascending=True)
        # Make round.inf with 4 decimals
        round_down_4_decimals = lambda x: int(x * 10 ** 4) / 10 ** 4
        output_df['round_down'] = output_df['lfn_biosample_replicate: N_ijk/N_jk'].apply(round_down_4_decimals)

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
            replicate_id_to_name = conn.execute(
                sqlalchemy.select([replicate_model.__table__.c.id, replicate_model.__table__.c.name])).fetchall()

        run_id_to_name_df = pandas.DataFrame.from_records(data=run_id_to_name, columns=['run_id', 'run_name'])
        output_df = output_df.merge(run_id_to_name_df, on='run_id')

        marker_id_to_name_df = pandas.DataFrame.from_records(data=marker_id_to_name, columns=['marker_id', 'marker_name'])
        output_df = output_df.merge(marker_id_to_name_df, on='marker_id')

        biosample_id_to_name_df = pandas.DataFrame.from_records(data=biosample_id_to_name, columns=['biosample_id', 'biosample_name'])
        output_df = output_df.merge(biosample_id_to_name_df, on='biosample_id')

        replicate_id_to_name_df = pandas.DataFrame.from_records(data=replicate_id_to_name, columns=['replicate_id', 'replicate_name'])
        output_df = output_df.merge(replicate_id_to_name_df, on='replicate_id')

        output_df = output_df[['run_name', 'marker_name', 'biosample_name', 'replicate_name', 'variant_id',
       'N_ijk', 'N_jk', 'lfn_biosample_replicate: N_ijk/N_jk', 'round_down']]

        ##########################################################
        #
        # 7. Write TSV file
        #
        ##########################################################

        output_df.sort_values(by=['lfn_biosample_replicate: N_ijk/N_jk', 'run_name', 'marker_name', 'biosample_name', 'replicate_name'],
                                       ascending=[True, True, True, True, True], inplace=True)
        output_df.to_csv(output_file_optimize_lfn, header=True, sep='\t', float_format='%.10f', index=False)

