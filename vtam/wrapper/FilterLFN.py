import sqlalchemy
from sqlalchemy import select
from vtam.utils.FilterLFNrunner import FilterLFNrunner
from vtam.utils.FastaInformation import FastaInformation
from vtam.utils.Logger import Logger
from vtam.utils.OptionManager import OptionManager
from vtam.utils.VTAMexception import VTAMexception
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper

import os
import pandas
import sys

from vtam.utils.VariantReadCountDF import VariantReadCountDF
from vtam.utils.VariantReadCountLikeTable import VariantReadCountLikeTable


class FilterLFN(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.FilterLFN"
    }

    # Input file
    __input_file_fastainfo = "fastainfo"
    # Input table
    __input_table_run = "Run"
    __input_table_marker = "Marker"
    __input_table_biosample = "Biosample"
    __input_table_replicate = "Replicate"
    __input_table_variant_read_count = "VariantReadCount"
    # Output table
    __output_table_filter_lfn = "FilterLFN"

    def specify_input_file(self):
        return[
            FilterLFN.__input_file_fastainfo,
        ]

    def specify_input_table(self):
        return [
            FilterLFN.__input_table_marker,
            FilterLFN.__input_table_run,
            FilterLFN.__input_table_biosample,
            FilterLFN.__input_table_replicate,
            FilterLFN.__input_table_variant_read_count,
        ]

    def specify_output_table(self):
        return [
            FilterLFN.__output_table_filter_lfn,
        ]

    def specify_params(self):
        return {
            "filter_lfn_variant": "required|int",
            "lfn_variant_threshold": "float",
            "lfn_variant_replicate_threshold": "float",
            "lfn_biosample_replicate_threshold": "required|float",
            "lfn_read_count_threshold": "required|float",
        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind

        ##########################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ##########################################################

        # Input file output
        input_file_fastainfo = self.input_file(FilterLFN.__input_file_fastainfo)
        # Add FilterLFNthresholdspecific
        # input_file_threshold_specific = self.input_file(FilterLFNthresholdspecific.__input_file_threshold_specific)
        #
        # Input table models
        run_model = self.input_table(FilterLFN.__input_table_run)
        marker_model = self.input_table(FilterLFN.__input_table_marker)
        biosample_model = self.input_table(FilterLFN.__input_table_biosample)
        replicate_model = self.input_table(FilterLFN.__input_table_replicate)
        input_variant_read_count_model = self.input_table(FilterLFN.__input_table_variant_read_count)
        #
        # Output table models
        output_filter_lfn_model = self.output_table(FilterLFN.__output_table_filter_lfn)
        #
        # Options
        filter_lfn_variant = int(self.option("filter_lfn_variant"))
        lfn_variant_threshold = self.option("lfn_variant_threshold")
        lfn_variant_replicate_threshold = self.option("lfn_variant_replicate_threshold")
        lfn_biosample_replicate_threshold = self.option("lfn_biosample_replicate_threshold")
        lfn_read_count_threshold = self.option("lfn_read_count_threshold")

        ##########################################################
        #
        # 1. Read fastainfo to get run_id, marker_id, biosample_id, replicate_id for current analysis
        # 2. Delete marker/run/biosample/replicate from variant_read_count_model
        # 3. Select marker/run/biosample/replicate from variant_read_count_model
        # 4. Apply filters
        # 5. Write filters to variant_filter_lfn_model
        #
        ##########################################################

        ##########################################################
        #
        # 1. Read fastainfo to get run_id, marker_id, biosample_id, replicate_id for current analysis
        #
        ##########################################################

        fasta_info = FastaInformation(input_file_fastainfo, engine, run_model, marker_model, biosample_model, replicate_model)
        fasta_info_record_list = fasta_info.get_fasta_information_record_list()

        ##########################################################
        #
        # 2. Delete marker/run/biosample/replicate from variant_read_count_model
        #
        ##########################################################

        variant_read_count_like_utils = VariantReadCountLikeTable(variant_read_count_like_model=output_filter_lfn_model, engine=engine)
        variant_read_count_like_utils.delete_output_filter_model(fasta_info_record_list=fasta_info_record_list)

        # with __engine.connect() as conn:
        #     conn.execute(filter_lfn_model.__table__.delete(), sample_instance_list)

        ##########################################################
        #
        #
        # 3. Select marker/run/biosample/replicate from variant_read_count_model
        #
        ##########################################################

        variant_read_count_df = fasta_info.get_variant_read_count_df(variant_read_count_like_model=input_variant_read_count_model, filter_id=None)

        ##########################################################
        #
        #
        ##########################################################

        #
        lfn_filter_runner = FilterLFNrunner(variant_read_count_df)
        #
        Logger.instance().info("Launching LFN filter:")
        #
        ############################################
        # Filter 2: f2_f4_lfn_delete_variant
        # Or
        # Filter  3: f3_f5_lfn_delete_variant_replicate
        ############################################
        if bool(filter_lfn_variant):
            lfn_filter_runner.f2_f4_lfn_delete_variant(lfn_variant_threshold)
        else:
            lfn_filter_runner.f3_f5_lfn_delete_variant_replicate(lfn_variant_replicate_threshold)
        #

        ############################################
        # Filter 6:  f6_lfn_delete_biosample_replicate_delete
        ############################################

        lfn_filter_runner.f6_lfn_delete_biosample_replicate(lfn_biosample_replicate_threshold)

        ############################################
        # Filter  7:f7_lfn_delete_absolute_read_count
        ############################################
        lfn_filter_runner.f7_lfn_delete_absolute_read_count(lfn_read_count_threshold)

        ############################################
        # Filter 8:f8_lfn_delete_do_not_pass_all_filters
        ############################################
        lfn_filter_runner.f8_lfn_delete_do_not_pass_all_filters()

        ############################################
        # Write to DB
        ############################################
        filter_output_df = lfn_filter_runner.variant_read_count_filter_delete_df

        record_list = VariantReadCountLikeTable.filter_delete_df_to_dict(filter_output_df)
        with engine.connect() as conn:
            try:
                conn.execute(output_filter_lfn_model.__table__.insert(), record_list)
            except sqlalchemy.exc.IntegrityError:
                Logger.instance().error(VTAMexception(
                    "Records are not unique"))

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


