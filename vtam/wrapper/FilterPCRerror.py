from wopmars.framework.database.tables.ToolWrapper import ToolWrapper

from vtam import Logger
from vtam.utils.FastaInformation import FastaInformation
from vtam.utils.FilterPCRerrorRunner import FilterPCRerrorRunner
from vtam.utils.VariantReadCountLikeTable import VariantReadCountLikeTable
from vtam.utils.OptionManager import OptionManager
from vtam.utils.PathManager import PathManager
from vtam.utils.VTAMexception import VTAMexception

import os
import pandas
import sys


class FilterPCRerror(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.FilterPCRerror"
    }

    # Input file
    __input_file_fastainfo = "fastainfo"
    # Input table
    __input_table_marker = "Marker"
    __input_table_run = "Run"
    __input_table_biosample = "Biosample"
    __input_table_replicate = "Replicate"
    __input_table_variant = "Variant"
    __input_table_filter_min_replicate_number = "FilterMinReplicateNumber"
    # Output table
    __output_table_filter_pcr_error = "FilterPCRerror"


    def specify_input_file(self):
        return[
            FilterPCRerror.__input_file_fastainfo,

        ]

    def specify_input_table(self):
        return [
            FilterPCRerror.__input_table_marker,
            FilterPCRerror.__input_table_run,
            FilterPCRerror.__input_table_biosample,
            FilterPCRerror.__input_table_replicate,
            FilterPCRerror.__input_table_variant,
            FilterPCRerror.__input_table_filter_min_replicate_number,
        ]


    def specify_output_table(self):
        return [
            FilterPCRerror.__output_table_filter_pcr_error,
        ]

    def specify_params(self):
        return {
            "pcr_error_var_prop": "float",
            # "log_verbosity": "int",
            # "log_file": "str"

        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        if not self.option("log_verbosity") is None:
            OptionManager.instance()['log_verbosity'] = int(self.option("log_verbosity"))
            OptionManager.instance()['log_file'] = str(self.option("log_file"))
        this_step_tmp_dir = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__))
        PathManager.mkdir_p(this_step_tmp_dir)

        ##########################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ##########################################################
        #
        # Input file output
        input_file_fastainfo = self.input_file(FilterPCRerror.__input_file_fastainfo)
        #
        # Input table models
        marker_model = self.input_table(FilterPCRerror.__input_table_marker)
        run_model = self.input_table(FilterPCRerror.__input_table_run)
        biosample_model = self.input_table(FilterPCRerror.__input_table_biosample)
        replicate_model = self.input_table(FilterPCRerror.__input_table_replicate)
        variant_model = self.input_table(FilterPCRerror.__input_table_variant)
        input_filter_min_replicate_model = self.input_table(FilterPCRerror.__input_table_filter_min_replicate_number)
        #
        # Options
        pcr_error_var_prop = self.option("pcr_error_var_prop")
        #
        # Output table models
        output_filter_pcr_error_model = self.output_table(FilterPCRerror.__output_table_filter_pcr_error)

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

        variant_read_count_like_utils = VariantReadCountLikeTable(variant_read_count_like_model=output_filter_pcr_error_model, engine=engine)
        variant_read_count_like_utils.delete_output_filter_model(fasta_info_record_list=fasta_info_record_list)

        ##########################################################
        #
        # variant_read_count_df
        #
        ##########################################################

        filter_id = None
        variant_read_count_df = fasta_info.get_variant_read_count_df(variant_read_count_like_model=input_filter_min_replicate_model, filter_id=filter_id)

        ##########################################################
        #
        # 4. Get variant_df: id, sequence
        #
        ##########################################################

        variant_df = fasta_info.get_variant_df(variant_read_count_like_model=input_filter_min_replicate_model,
                                                          variant_model=variant_model)

        ################################################################################################################
        #
        # Run per biosample_id
        #
        ################################################################################################################

        record_list = []

        run_marker_biosample_df = variant_read_count_df[['run_id', 'marker_id', 'biosample_id']].drop_duplicates()
        for row in run_marker_biosample_df.itertuples():
            run_id = row.run_id
            marker_id = row.marker_id
            biosample_id = row.biosample_id
            # for biosample_id in variant_read_count_df.biosample_id.unique().tolist():

            # variant_read_count_per_biosample_df = variant_read_count_df.loc[variant_read_count_df.biosample_id==biosample_id]
            variant_read_count_per_biosample_df = variant_read_count_df.loc[(variant_read_count_df.run_id == run_id)
                                                                   & (variant_read_count_df.marker_id == marker_id)
                                                                   & (variant_read_count_df.biosample_id == biosample_id)]

            variant_per_biosample_df = variant_df.loc[variant_df.id.isin(variant_read_count_df.variant_id.unique().tolist())]
            this_step_tmp_per_biosample_dir = os.path.join(this_step_tmp_dir, "run_{}_marker_{}_biosample{}".format(run_id, marker_id, biosample_id))
            PathManager.instance().mkdir_p(this_step_tmp_per_biosample_dir)

            ##########################################################
            #
            # Run vsearch and get alignement df
            #
            ##########################################################

            filter_pcr_error_runner = FilterPCRerrorRunner(variant_expected_df=variant_per_biosample_df,
                                                           variant_unexpected_df=variant_per_biosample_df,
                                                           variant_read_count_df=variant_read_count_per_biosample_df,
                                                           tmp_dir=this_step_tmp_per_biosample_dir)
            filter_output_per_biosample_df = filter_pcr_error_runner.get_filter_output_df(pcr_error_var_prop)

            ##########################################################
            #
            # Per biosample add to record list
            #
            ##########################################################

            record_per_biosample_list = VariantReadCountLikeTable.filter_delete_df_to_dict(filter_output_per_biosample_df)
            record_list = record_list + record_per_biosample_list

        ##########################################################
        #
        # Write to DB
        #
        ##########################################################

        with engine.connect() as conn:
            conn.execute(output_filter_pcr_error_model.__table__.insert(), record_list)

        ##########################################################
        #
        # Exit vtam if all variants delete
        #
        ##########################################################

        filter_output_df = pandas.DataFrame.from_records(data=record_list)

        try:
            assert not filter_output_df.filter_delete.sum() == filter_output_df.shape[0]
        except AssertionError:
            Logger.instance().warning(VTAMexception("This filter has deleted all the variants: {}. The analysis will stop here.".format(self.__class__.__name__)))
            sys.exit(0)

