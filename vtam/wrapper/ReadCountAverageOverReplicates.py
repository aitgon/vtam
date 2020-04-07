from sqlalchemy import select, bindparam

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
    __input_file_readinfo = "readinfo"
    __input_table_filter_codon_stop = "FilterCodonStop"
    # Output table
    __output_table_filter_consensus = "ReadCountAverageOverReplicates"



    def specify_input_file(self):
        return[
            ReadCountAverageOverReplicates.__input_file_readinfo,

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
        input_file_readinfo = self.input_file(ReadCountAverageOverReplicates.__input_file_readinfo)
        #
        # Input table models
        marker_model = self.input_table(ReadCountAverageOverReplicates.__input_table_marker)
        run_model = self.input_table(ReadCountAverageOverReplicates.__input_table_run)
        codon_stop_model = self.input_table(ReadCountAverageOverReplicates.__input_table_filter_codon_stop)
        biosample_model = self.input_table(ReadCountAverageOverReplicates.__input_table_biosample)

        #
        # Output table models
        consensus_model = self.output_table(ReadCountAverageOverReplicates.__output_table_filter_consensus)

        ##########################################################
        #
        # 1. Read readinfo to get run_id, marker_id, biosample_id, replicate for current analysis
        #
        ##########################################################

        fasta_info_tsv = FastaInformationTSV(engine=engine, fasta_info_tsv=input_file_readinfo)

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


        ##########################################################
        #
        # 3. Select marker/run/biosample/replicate from variant_read_count_model
        #
        ##########################################################

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

        ############################################
        # Write to DB
        ############################################
        record_list = VariantReadCountLikeTable.filter_delete_df_to_dict(df_out)
        with engine.connect() as conn:

            # Delete instances that will be inserted
            del_stmt = consensus_model.__table__.delete() \
                .where(consensus_model.run_id == bindparam('run_id')) \
                .where(consensus_model.marker_id == bindparam('marker_id')) \
                .where(consensus_model.biosample_id == bindparam('biosample_id')) \
                .where(consensus_model.replicate == bindparam('replicate'))
            conn.execute(del_stmt, record_list)

            # Insert new instances
            conn.execute(consensus_model.__table__.insert(), record_list)

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
