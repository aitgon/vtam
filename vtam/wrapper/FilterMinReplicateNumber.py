from vtam.utils.RunnerFilterMinReplicateNumber import RunnerFilterMinReplicateNumber
from vtam.utils.FileSampleInformation import FileSampleInformation
from vtam.utils.DataframeVariantReadCountLike import DataframeVariantReadCountLike
from vtam.utils.ModelVariantReadCountLike import ModelVariantReadCountLike
from vtam.utils.Logger import Logger
from vtam.utils.VTAMexception import VTAMexception
from wopmars.models.ToolWrapper import ToolWrapper

import pandas
import sys


class FilterMinReplicateNumber(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.FilterMinReplicateNumber"
    }

    # Input file
    __input_file_sortedinfo = "sortedinfo"
    # Input table
    __input_table_run = "Run"
    __input_table_marker = "Marker"
    __input_table_sample = "Sample"
    __input_table_variant_filter_lfn = "FilterLFN"
    # Output table
    __output_table_filter_min_replicate_number = "FilterMinReplicateNumber"

    def specify_input_file(self):
        return[
            "sortedinfo", "params",
        ]

    def specify_input_table(self):
        return [
            FilterMinReplicateNumber.__input_table_run,
            FilterMinReplicateNumber.__input_table_marker,
            FilterMinReplicateNumber.__input_table_sample,
            FilterMinReplicateNumber.__input_table_variant_filter_lfn,
        ]

    def specify_output_table(self):
        return [
            FilterMinReplicateNumber.__output_table_filter_min_replicate_number, ]

    def specify_params(self):
        return {
            "min_replicate_number": "int",
        }

    def run(self):
        session = self.session
        engine = session._session().get_bind()

        #######################################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        #######################################################################
        #
        # Input files
        fasta_info_tsv = self.input_file(
            FilterMinReplicateNumber.__input_file_sortedinfo)
        #
        # Input tables
        input_filter_lfn_model = self.input_table(
            FilterMinReplicateNumber.__input_table_variant_filter_lfn)
        #
        # Options
        min_replicate_number = self.option("min_replicate_number")
        # input_filter_lfn = self.option("input_filter_lfn")
        #
        # Output tables
        output_filter_min_replicate_model = self.output_table(
            FilterMinReplicateNumber.__output_table_filter_min_replicate_number)

        #######################################################################
        #
        # 1. Read sortedinfo to get run_id, marker_id, sample_id, replicate for current analysis
        # 2. Delete marker_name/run_name/sample/replicate from variant_read_count_model
        # 3. Get nijk_df input
        #
        #######################################################################

        sample_info_tsv_obj = FileSampleInformation(tsv_path=fasta_info_tsv)

        sample_info_tsv_obj.delete_from_db(
            engine=engine, variant_read_count_like_model=output_filter_min_replicate_model)
        filter_id = None
        if input_filter_lfn_model.__tablename__ == "FilterLFN":
            filter_id = 8  # Variant pass all filters LFN
        variant_read_count_df = sample_info_tsv_obj.get_nijk_df(
            variant_read_count_like_model=input_filter_lfn_model, engine=engine, filter_id=filter_id)

        #######################################################################
        #
        # 4. Run Filter
        #
        #######################################################################

        variant_read_count_delete_df = RunnerFilterMinReplicateNumber(
            variant_read_count_df) .get_variant_read_count_delete_df(min_replicate_number)

        #######################################################################
        #
        # 5. Write to DB
        # 6. Touch output tables, to update modification date
        # 7. Exit vtam if all variants delete
        #
        #######################################################################

        DataframeVariantReadCountLike(variant_read_count_delete_df).to_sql(
            engine=engine, variant_read_count_like_model=output_filter_min_replicate_model)

        for output_table_i in self.specify_output_table():
            declarative_meta_i = self.output_table(output_table_i)
            obj = session.query(declarative_meta_i).order_by(
                declarative_meta_i.id.desc()).first()
            session.query(declarative_meta_i).filter_by(
                id=obj.id).update({'id': obj.id})
            session.commit()

        if variant_read_count_delete_df.filter_delete.sum(
        ) == variant_read_count_delete_df.shape[0]:
            Logger.instance().warning(
                VTAMexception(
                    "This filter has deleted all the variants: {}. "
                    "The analysis will stop here.".format(
                        self.__class__.__name__)))
            sys.exit(0)
