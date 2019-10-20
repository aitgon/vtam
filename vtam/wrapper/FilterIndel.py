from sqlalchemy import select
from vtam.utils.FilterCommon import FilterCommon
from vtam.utils.Logger import Logger
from vtam.utils.OptionManager import OptionManager
from vtam.utils.VTAMexception import VTAMexception
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper

import pandas
import sys


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
        biosample_model = self.input_table(FilterIndel.__input_table_biosample)
        replicate_model = self.input_table(FilterIndel.__input_table_replicate)
        variant_model = self.input_table(FilterIndel.__input_table_Variant)
        input_filter_model = self.input_table(FilterIndel.__input_table_filter_renkonen)
        #
        # Output table models
        output_filter_models = self.output_table(FilterIndel.__output_table_filter_indel)

        ##########################################################
        #
        # 1. Read fastainfo to get run_id, marker_id, biosample_id, replicate_id for current analysis
        #
        ##########################################################

        filter_various = FilterCommon(self.__class__.__name__, engine, run_model, marker_model, biosample_model, replicate_model,
                                      input_filter_model,
                                      output_filter_models=output_filter_models)
        fastainfo_instance_list = filter_various.get_fastainfo_instance_list_with_ids(input_file_fastainfo)

        ##########################################################
        #
        # 2. Delete /run/markerbiosample/replicate from this filter table
        #
        ##########################################################

        filter_various.delete_output_filter_model(fastainfo_instance_list)


        ##########################################################
        #
        # 3. Select variant_read_count_model
        #
        ##########################################################

        variant_read_count_df = filter_various.get_variant_read_count_model(fastainfo_instance_list)

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

        ##########################################################
        #
        # Write to DB
        #
        ##########################################################

        records = FilterCommon.filter_delete_df_to_dict(filter_output_df)
        with engine.connect() as conn:
            conn.execute(output_filter_models.__table__.insert(), records)

        ##########################################################
        #
        # Exit vtam if all variants delete
        #
        ##########################################################

        try:
            assert not filter_output_df.filter_delete.sum() == filter_output_df.shape[0]
        except AssertionError:
            Logger.instance().warning(VTAMexception("This filter has deleted all the variants: {}. The analysis will stop here.".format(self.__class__.__name__)))
            sys.exit(0)




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


