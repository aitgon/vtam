from sqlalchemy import select

from vtam.utils.SampleInformationUtils import FastaInformationTSV
from vtam.utils.VariantReadCountLikeTable import VariantReadCountLikeTable
from wopmars.models.ToolWrapper import ToolWrapper

import os
import pandas
import sys


class ReadCountAverageOverReplicates(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.ReadCountAverageOverReplicates"
    }

    # Input file
    # Input table
    __input_table_marker = "Marker"
    __input_table_run = "Run"
    __input_table_biosample = "Biosample"
    __input_file_fastainfo = "fastainfo"
    __input_table_filter_codon_stop = "FilterCodonStop"
    # Output table
    __output_table_filter_consensus = "ReadCountAverageOverReplicates"



    def specify_input_file(self):
        return[
            ReadCountAverageOverReplicates.__input_file_fastainfo,

        ]

    def specify_input_table(self):
        return [
            ReadCountAverageOverReplicates.__input_table_marker,
            ReadCountAverageOverReplicates.__input_table_run,
            ReadCountAverageOverReplicates.__input_table_biosample,
            ReadCountAverageOverReplicates.__input_table_filter_codon_stop,

        ]


    def specify_output_table(self):
        return [
            ReadCountAverageOverReplicates.__output_table_filter_consensus,

        ]

    def specify_params(self):
        return {
        }



    def run(self):
        session = self.session
        engine = session._session().get_bind()
        #
        # Input file output
        input_file_fastainfo = self.input_file(ReadCountAverageOverReplicates.__input_file_fastainfo)
        #
        # Input table models
        marker_model = self.input_table(ReadCountAverageOverReplicates.__input_table_marker)
        run_model = self.input_table(ReadCountAverageOverReplicates.__input_table_run)
        codon_stop_model = self.input_table(ReadCountAverageOverReplicates.__input_table_filter_codon_stop)
        biosample_model = self.input_table(ReadCountAverageOverReplicates.__input_table_biosample)

        #
        # Output table models
        consensus_model = self.output_table(ReadCountAverageOverReplicates.__output_table_filter_consensus)


        # ##########################################################
        # #
        # # 1. Read fastainfo to get run_id, marker_id, biosample_id, replicate for current analysis
        # #
        # ##########################################################
        # fastainfo_df = pandas.read_csv(input_file_fastainfo, sep="\t", header=0,\
        #     names=['tag_forward', 'primer_forward', 'tag_reverse', 'primer_reverse', 'marker_name', 'biosample_name',\
        #     'replicate', 'run_name', 'fastq_fwd', 'fastq_rev', 'fasta_path'])
        # sample_instance_list = []
        # for row in fastainfo_df.itertuples():
        #     marker_name = row.marker_name
        #     run_name = row.run_name
        #     biosample_name = row.biosample_name
        #     replicate = row.replicate
        #     with engine.connect() as conn:
        #         # get run_id ###########
        #         stmt_select_run_id = select([run_model.__table__.c.id]).where(run_model.__table__.c.name==run_name)
        #         run_id = conn.execute(stmt_select_run_id).first()[0]
        #         # get marker_id ###########
        #         stmt_select_marker_id = select([marker_model.__table__.c.id]).where(marker_model.__table__.c.name==marker_name)
        #         marker_id = conn.execute(stmt_select_marker_id).first()[0]
        #         # get biosample_id ###########
        #         stmt_select_biosample_id = select([biosample_model.__table__.c.id]).where(biosample_model.__table__.c.name==biosample_name)
        #         biosample_id = conn.execute(stmt_select_biosample_id).first()[0]
        #         # add this sample_instance ###########
        #         sample_instance_list.append({'run_id': run_id, 'marker_id': marker_id, 'biosample_id':biosample_id, 'replicate':replicate})

        ##########################################################
        #
        # 1. Read fastainfo to get run_id, marker_id, biosample_id, replicate for current analysis
        #
        ##########################################################

        fasta_info_tsv = FastaInformationTSV(engine=engine, fasta_info_tsv=input_file_fastainfo, run_model=run_model,
                                             marker_model=marker_model, biosample_model=biosample_model)

        ##########################################################
        #
        # 2. Delete /run/markerbiosample/replicate from this filter table
        #
        ##########################################################
        # with engine.connect() as conn:
        #     # conn.execute(consensus_model.__table__.delete(), sample_instance_list)
        #     conn.execute(consensus_model.__table__.delete(), sample_instance_list)
        #
        variant_read_count_like_utils = VariantReadCountLikeTable(
            variant_read_count_like_model=consensus_model, engine=engine)
        variant_read_count_like_utils.delete_from_db(sample_record_list=fasta_info_tsv.sample_record_list)


        # ##########################################################
        # #
        # # 3. Select marker/run/biosample/replicate from variant_read_count_model
        # #
        # ##########################################################
        #
        # codon_stop_model_table = codon_stop_model.__table__
        #
        # variant_read_count_list = []
        # for sample_instance in sample_instance_list:
        #     run_id = sample_instance['run_id']
        #     marker_id = sample_instance['marker_id']
        #     biosample_id = sample_instance['biosample_id']
        #     replicate = sample_instance['replicate']
        #     stmt_select = select([codon_stop_model_table.c.run_id,
        #                           codon_stop_model_table.c.marker_id,
        #                           codon_stop_model_table.c.biosample_id,
        #                           codon_stop_model_table.c.replicate,
        #                           codon_stop_model_table.c.variant_id,
        #                           codon_stop_model_table.c.read_count]).distinct()\
        #                             .where(codon_stop_model_table.c.run_id == run_id)\
        #                             .where(codon_stop_model_table.c.marker_id == marker_id)\
        #                             .where(codon_stop_model_table.c.biosample_id == biosample_id)\
        #                             .where(codon_stop_model_table.c.replicate == replicate)\
        #                             .where(codon_stop_model_table.c.filter_delete == 0)
        #     with engine.connect() as conn:
        #         for row2 in conn.execute(stmt_select).fetchall():
        #             variant_read_count_list.append(row2)
        # #
        # variant_read_count_df = pandas.DataFrame.from_records(variant_read_count_list,
        #     columns=['run_id', 'marker_id', 'biosample_id', 'replicate', 'variant_id', 'read_count'])

        variant_read_count_df = fasta_info_tsv.get_variant_read_count_df(
            variant_read_count_like_model=codon_stop_model, filter_id=None)

        # Exit if no variants for analysis
        try:
            assert variant_read_count_df.shape[0] > 0
        except AssertionError:
            sys.stderr.write("Error: No variants available for this filter: {}".format(os.path.basename(__file__)))
            sys.exit(1)

        ##########################################################
        #
        #
        ##########################################################
        # else:
        ##########################################################
        #
        # 4. Run Filter
        #
        ##########################################################
        df_out = read_count_average_over_replicates(variant_read_count_df)

        # ##########################################################
        # #
        # # 5. Insert Filter data
        # #
        # ##########################################################
        # records = df_out.to_dict('records')
        # with __engine.connect() as conn:
        #         conn.execute(consensus_model.__table__.insert(), records)

        ############################################
        # Write to DB
        ############################################
        records = VariantReadCountLikeTable.filter_delete_df_to_dict(df_out)
        with engine.connect() as conn:
            conn.execute(consensus_model.__table__.insert(), records)



def read_count_average_over_replicates(variant_read_count_df):
    """
        Function used to display the read average of the remaining variant
    """
    variants_sequences = variant_read_count_df['variant_id'].copy()
    variants_sequences = list(set(variants_sequences))
    read_average_columns = ['variant', 'read_average']
    read_average_df = pandas.DataFrame(columns=read_average_columns)

    # sum of read_count over variant_id and biosample_id
    read_count_sum_over_variant_id_and_biosample_id_df = variant_read_count_df.groupby(['run_id', 'marker_id', 'variant_id', 'biosample_id']).sum().reset_index()
    read_count_sum_over_variant_id_and_biosample_id_df.drop('replicate', axis=1, inplace=True)
    read_count_sum_over_variant_id_and_biosample_id_df = read_count_sum_over_variant_id_and_biosample_id_df.rename(columns={'read_count': 'read_count'})

    #  count of replicate number per variant_id and biosample_id
    replicate_count_over_variant_id_and_biosample_id_df = variant_read_count_df.groupby(['run_id', 'marker_id', 'variant_id', 'biosample_id']).count().reset_index()
    replicate_count_over_variant_id_and_biosample_id_df.drop('read_count', axis=1, inplace=True)
    replicate_count_over_variant_id_and_biosample_id_df = replicate_count_over_variant_id_and_biosample_id_df.rename(columns={'replicate': 'replicate_count'})

    # merge
    df_out = read_count_sum_over_variant_id_and_biosample_id_df.merge(replicate_count_over_variant_id_and_biosample_id_df, left_on=('run_id', 'marker_id', 'variant_id', 'biosample_id'),right_on=('run_id', 'marker_id', 'variant_id', 'biosample_id'))
    df_out['read_count_average'] = df_out.read_count/df_out.replicate_count
    #
    return df_out



