import inspect
import os

import sys
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper


from sqlalchemy import select
import pandas

from vtam.utils.OptionManager import OptionManager
from vtam.utils.VTAMexception import VTAMexception
from vtam.utils.Logger import Logger
from vtam.utils.utilities import filter_delete_df_to_dict


class FilterIndel(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.FilterIndel"
    }

    # Input file
    __input_file_fastainfo = "fastainfo"
    # Input table
    __input_table_marker = "Marker"
    __input_table_run = "Run"
    __input_table_biosample = "Biosample"
    __input_table_replicate = "Replicate"
    __input_table_filter_renkonen = "FilterRenkonen"
    __input_table_Variant = "Variant"
    # Output table
    __output_table_filter_indel = "FilterIndel"



    def specify_input_file(self):
        return[
            FilterIndel.__input_file_fastainfo,

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

        ##########################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ##########################################################
        #
        # Input file output
        input_file_fastainfo = self.input_file(FilterIndel.__input_file_fastainfo)
        #
        # Input table models
        marker_model = self.input_table(FilterIndel.__input_table_marker)
        run_model = self.input_table(FilterIndel.__input_table_run)
        renkonen_model = self.input_table(FilterIndel.__input_table_filter_renkonen)
        biosample_model = self.input_table(FilterIndel.__input_table_biosample)
        replicate_model = self.input_table(FilterIndel.__input_table_replicate)
        variant_model = self.input_table(FilterIndel.__input_table_Variant)
        #
        # Output table models
        indel_model = self.output_table(FilterIndel.__output_table_filter_indel)

        ##########################################################
        #
        # 2. Read fastainfo to get run_id, marker_id, biosample_id, replicate_id for current analysis
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
        # 3. Delete /run/markerbiosample/replicate from this filter table
        #
        ##########################################################
        with engine.connect() as conn:
            conn.execute(indel_model.__table__.delete(), sample_instance_list)
        #


        ##########################################################
        #
        # 4. Select marker/run/biosample/replicate from variant_read_count_model
        #
        ##########################################################

        renkonen_model_table = renkonen_model.__table__

        variant_read_count_list = []
        for sample_instance in sample_instance_list:
            run_id = sample_instance['run_id']
            marker_id = sample_instance['marker_id']
            biosample_id = sample_instance['biosample_id']
            replicate_id = sample_instance['replicate_id']
            stmt_select = select([renkonen_model_table.c.run_id,
                                  renkonen_model_table.c.marker_id,
                                  renkonen_model_table.c.biosample_id,
                                  renkonen_model_table.c.replicate_id,
                                  renkonen_model_table.c.variant_id,
                                  renkonen_model_table.c.read_count]).distinct()\
                                    .where(renkonen_model_table.c.run_id == run_id)\
                                    .where(renkonen_model_table.c.marker_id == marker_id)\
                                    .where(renkonen_model_table.c.biosample_id == biosample_id)\
                                    .where(renkonen_model_table.c.replicate_id == replicate_id)\
                                    .where(renkonen_model_table.c.filter_delete == 0)
            with engine.connect() as conn:
                for row2 in conn.execute(stmt_select).fetchall():
                    variant_read_count_list.append(row2)
        #
        variant_read_count_df = pandas.DataFrame.from_records(variant_read_count_list,
            columns=['run_id', 'marker_id', 'biosample_id', 'replicate_id', 'variant_id', 'read_count'])

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
        # run_id, marker_id, variant_id, biosample_id, replicate_id, read_count, filter_delete
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
        # 5. Run Filter
        #
        ##########################################################
        filter_output_df = f13_filter_indel(variant_read_count_df, variant_df)

        # ##########################################################
        # #
        # # 6. Insert Filter data
        # #
        # ##########################################################
        #
        # records = df_out.to_dict('records')
        # with engine.connect() as conn:
        #         conn.execute(indel_model.__table__.insert(), records)

        ############################################
        # Write to DB
        ############################################
        records = filter_delete_df_to_dict(filter_output_df)
        with engine.connect() as conn:
            conn.execute(indel_model.__table__.insert(), records)

        ##########################################################
        #
        # 6. Exit vtam if all variants delete
        #
        ##########################################################
        # Exit if no variants for analysis
        try:
            assert not filter_output_df.shape[0] == 0
        except AssertionError:
            Logger.instance().info(VTAMexception("Error: This filter has deleted all the variants"))
            sys.exit(1)




def f13_filter_indel(variant_read_count_df, variant_df):
    """
    filter chimera
    """

    df_out = variant_read_count_df.copy()
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


