import inspect

from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmetabarcoding.utils.Logger import Logger
from wopmetabarcoding.utils.OptionManager import OptionManager


from sqlalchemy import select
import pandas



class ReadCountAverageOverReplicates(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.ReadCountAverageOverReplicates"
    }

    # Input file
    # Input table
    __input_table_marker = "Marker"
    __input_table_run = "Run"
    __input_table_biosample = "Biosample"
    __input_table_replicate = "Replicate"
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
            ReadCountAverageOverReplicates.__input_table_replicate,
            ReadCountAverageOverReplicates.__input_table_filter_codon_stop,

        ]


    def specify_output_table(self):
        return [
            ReadCountAverageOverReplicates.__output_table_filter_consensus,

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
        #
        # Input file path
        input_file_fastainfo = self.input_file(ReadCountAverageOverReplicates.__input_file_fastainfo)
        #
        # Input table models
        marker_model = self.input_table(ReadCountAverageOverReplicates.__input_table_marker)
        run_model = self.input_table(ReadCountAverageOverReplicates.__input_table_run)
        codon_stop_model = self.input_table(ReadCountAverageOverReplicates.__input_table_filter_codon_stop)
        biosample_model = self.input_table(ReadCountAverageOverReplicates.__input_table_biosample)
        replicate_model = self.input_table(ReadCountAverageOverReplicates.__input_table_replicate)

        #

        #
        # Output table models
        consensus_model = self.output_table(ReadCountAverageOverReplicates.__output_table_filter_consensus)


        ##########################################################
        #
        # 1. Read fastainfo to get run_id, marker_id, biosample_id, replicate_id for current analysis
        #
        ##########################################################
        fastainfo_df = pandas.read_csv(input_file_fastainfo, sep="\t", header=0,\
            names=['tag_forward', 'primer_forward', 'tag_reverse', 'primer_reverse', 'marker_name', 'biosample_name',\
            'replicate_name', 'run_name', 'fastq_fwd', 'fastq_rev', 'fasta'])
        sample_instance_list = []
        for row in fastainfo_df.itertuples():
            marker_name = row.marker_name
            run_name = row.run_name
            biosample_name = row.biosample_name
            replicate_name = row.replicate_name
            with engine.connect() as conn:
                # get run_id ###########
                stmt_select_run_id = select([run_model.__table__.c.id]).where(run_model.__table__.c.name==run_name)
                run_id = conn.execute(stmt_select_run_id).first()[0]
                # get marker_id ###########
                stmt_select_marker_id = select([marker_model.__table__.c.id]).where(marker_model.__table__.c.name==marker_name)
                marker_id = conn.execute(stmt_select_marker_id).first()[0]
                # get biosample_id ###########
                stmt_select_biosample_id = select([biosample_model.__table__.c.id]).where(biosample_model.__table__.c.name==biosample_name)
                biosample_id = conn.execute(stmt_select_biosample_id).first()[0]
                # get replicate_id ###########
                stmt_select_replicate_id = select([replicate_model.__table__.c.id]).where(replicate_model.__table__.c.name==replicate_name)
                replicate_id = conn.execute(stmt_select_replicate_id).first()[0]
                # add this sample_instance ###########
                sample_instance_list.append({'run_id': run_id, 'marker_id': marker_id, 'biosample_id':biosample_id, 'replicate_id':replicate_id})

        ##########################################################
        #
        # 2. Delete /run/markerbiosample/replicate from this filter table
        #
        ##########################################################
        with engine.connect() as conn:
            # conn.execute(consensus_model.__table__.delete(), sample_instance_list)
            conn.execute(consensus_model.__table__.delete(), sample_instance_list)
        #


        ##########################################################
        #
        # 3. Select marker/run/biosample/replicate from variant_read_count_model
        #
        ##########################################################

        codon_stop_model_table = codon_stop_model.__table__
        stmt_variant_filter_lfn = select([codon_stop_model_table.c.marker_id,
                                          codon_stop_model_table.c.run_id,
                                          codon_stop_model_table.c.variant_id,
                                          codon_stop_model_table.c.biosample_id,
                                          codon_stop_model_table.c.replicate_id,
                                          codon_stop_model_table.c.read_count])\
            .where(codon_stop_model_table.c.filter_delete == 0)

        # Select to DataFrame
        variant_filter_lfn_passed_list = []
        with engine.connect() as conn:
            for row in conn.execute(stmt_variant_filter_lfn).fetchall():
                variant_filter_lfn_passed_list.append(row)
        variant_read_count_df = pandas.DataFrame.from_records(variant_filter_lfn_passed_list,
                    columns=['marker_id','run_id', 'variant_id', 'biosample_id', 'replicate_id', 'read_count'])
        if variant_read_count_df.shape[0] == 0:
            Logger.instance().debug(
                "file: {}; line: {}; No data input for this filter.".format(__file__,
                                                                      inspect.currentframe().f_lineno,
                                                                      'Consensus'))
        else:
            ##########################################################
            #
            # 4. Run Filter
            #
            ##########################################################
            df_out = read_count_average_over_replicates(variant_read_count_df)

            ##########################################################
            #
            # 5. Insert Filter data
            #
            ##########################################################
            records = df_out.to_dict('records')
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
    read_count_sum_over_variant_id_and_biosample_id_df.drop('replicate_id', axis=1, inplace=True)
    read_count_sum_over_variant_id_and_biosample_id_df = read_count_sum_over_variant_id_and_biosample_id_df.rename(columns={'read_count': 'read_count'})

    #  count of replicate number per variant_id and biosample_id
    replicate_count_over_variant_id_and_biosample_id_df = variant_read_count_df.groupby(['run_id', 'marker_id', 'variant_id', 'biosample_id']).count().reset_index()
    replicate_count_over_variant_id_and_biosample_id_df.drop('read_count', axis=1, inplace=True)
    replicate_count_over_variant_id_and_biosample_id_df = replicate_count_over_variant_id_and_biosample_id_df.rename(columns={'replicate_id': 'replicate_count'})

    # merge
    df_out = read_count_sum_over_variant_id_and_biosample_id_df.merge(replicate_count_over_variant_id_and_biosample_id_df, left_on=('run_id', 'marker_id', 'variant_id', 'biosample_id'),right_on=('run_id', 'marker_id', 'variant_id', 'biosample_id'))
    df_out['read_count_average'] = df_out.read_count/df_out.replicate_count
    #
    return df_out



