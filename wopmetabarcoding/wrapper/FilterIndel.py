import inspect

import sys
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper


from sqlalchemy import select
import pandas

from wopmetabarcoding.utils.logger import logger


class FilterIndel(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.FilterIndel"
    }

    # Input file
    # Input table
    __input_table_marker = "Marker"
    __input_table_run = "Run"
    __input_table_biosample = "Biosample"
    __input_table_replicate = "Replicate"
    __input_file_sample2fasta = "sample2fasta"
    __input_table_filter_renkonen = "FilterRenkonen"
    __input_table_Variant = "Variant"
    # Output table
    __output_table_filter_indel = "FilterIndel"



    def specify_input_file(self):
        return[
            FilterIndel.__input_file_sample2fasta,

        ]

    def specify_input_table(self):
        return [
            FilterIndel.__input_table_marker,
            FilterIndel.__input_table_run,
            FilterIndel.__input_table_biosample,
            FilterIndel.__input_table_replicate,
            FilterIndel.__input_table_filter_renkonen,
            FilterIndel.__input_table_Variant,
        ]


    def specify_output_table(self):
        return [
            FilterIndel.__output_table_filter_indel,

        ]



    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        #
        # Input file path
        input_file_sample2fasta = self.input_file(FilterIndel.__input_file_sample2fasta)
        #
        # Input table models
        marker_model = self.input_table(FilterIndel.__input_table_marker)
        run_model = self.input_table(FilterIndel.__input_table_run)
        renkonen_model = self.input_table(FilterIndel.__input_table_filter_renkonen)
        biosample_model = self.input_table(FilterIndel.__input_table_biosample)
        replicate_model = self.input_table(FilterIndel.__input_table_replicate)
        variant_model = self.input_table(FilterIndel.__input_table_Variant)
        #

        #
        # Output table models
        indel_model = self.output_table(FilterIndel.__output_table_filter_indel)


        ##########################################################
        #
        # 1. Read sample2fasta to get run_id, marker_id, biosample_id, replicate_id for current analysis
        #
        ##########################################################
        sample2fasta_df = pandas.read_csv(input_file_sample2fasta, sep="\t", header=None,\
            names=['tag_forward', 'primer_forward', 'tag_reverse', 'primer_reverse', 'marker_name', 'biosample_name',\
            'replicate_name', 'run_name', 'fastq_fwd', 'fastq_rev', 'fasta'])
        sample_instance_list = []
        for row in sample2fasta_df.itertuples():
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
            conn.execute(indel_model.__table__.delete(), sample_instance_list)
        #


        ##########################################################
        #
        # 3. Select marker/run/biosample/replicate from variant_read_count_model
        #
        ##########################################################

        renkonen_model_table = renkonen_model.__table__
        stmt_variant_filter_lfn = select([renkonen_model_table.c.marker_id,
                                          renkonen_model_table.c.run_id,
                                          renkonen_model_table.c.variant_id,
                                          renkonen_model_table.c.biosample_id,
                                          renkonen_model_table.c.replicate_id,
                                          # renkonen_model_table.c.filter_id,
                                          # renkonen_model_table.c.filter_delete,
                                          renkonen_model_table.c.read_count])\
            .where(renkonen_model_table.c.filter_id == 12)\
            .where(renkonen_model_table.c.filter_delete == 0)
        # Select to DataFrame
        variant_filter_lfn_passed_list = []
        with engine.connect() as conn:
            for row in conn.execute(stmt_variant_filter_lfn).fetchall():
                variant_filter_lfn_passed_list.append(row)
        variant_read_count_df = pandas.DataFrame.from_records(variant_filter_lfn_passed_list,
                    columns=['marker_id','run_id', 'variant_id', 'biosample_id', 'replicate_id', 'read_count'])
        if variant_read_count_df.shape[0] != 0:
            logger.debug(
                "file: {}; line: {}; No data input for this filter.".format(__file__,
                                                                      inspect.currentframe().f_lineno,
                                                                      'Indel'))

            # run_id, marker_id, variant_id, biosample_id, replicate_id, read_count, filter_id, filter_delete
            variant_model_table = variant_model.__table__
            stmt_variant = select([variant_model_table.c.id,
                                   variant_model_table.c.sequence])

            # Select to DataFrame
            variant_filter_lfn_passed_list = []
            with engine.connect() as conn:
                for row in conn.execute(stmt_variant).fetchall():
                    variant_filter_lfn_passed_list.append(row)
            variant_df = pandas.DataFrame.from_records(variant_filter_lfn_passed_list,
                                                                  columns=['id', 'sequence'])
            ##########################################################
            #
            # 4. Run Filter
            #
            ##########################################################
            df_out = f13_filter_indel(variant_read_count_df, variant_df)

            ##########################################################
            #
            # 5. Insert Filter data
            #
            ##########################################################

            records = df_out.to_dict('records')
            with engine.connect() as conn:
                    conn.execute(indel_model.__table__.insert(), records)




def f13_filter_indel(variant_read_count_df, variant_df):
    """
    filter chimera
    """

    df_out = variant_read_count_df.copy()
    df_out['filter_id'] = 13
    df_out['filter_delete'] = False
    #
    df = variant_df.copy()
    df['sequence_length_module_3'] = variant_df.sequence.apply(lambda x: len(x) % 3) # compute module for each variant
    majority_sequence_length_module_3 = df.sequence_length_module_3.mode() # most common remaining of modulo 3
    # select id of variant that do not pass on a list
    df = df.loc[df['sequence_length_module_3'] != majority_sequence_length_module_3.values[0]]
    do_not_pass_variant_id_list = df.id.tolist()
    #
    for id in do_not_pass_variant_id_list:
        df_out.loc[df_out['variant_id'] == id, 'filter_delete'] = True
    #
    return df_out


