import pathlib

from vtam.utils.FileCutoffSpecific import FileCutoffSpecific

from vtam.utils.RunnerFilterLFN import RunnerFilterLFN
from vtam.utils.Logger import Logger
from vtam.utils.FileSampleInformation import FileSampleInformation
from vtam.utils.VTAMexception import VTAMexception
from wopmars.models.ToolWrapper import ToolWrapper

import sys

from vtam.utils.DataframeVariantReadCountLike import DataframeVariantReadCountLike
from vtam.utils.ModelVariantReadCountLike import ModelVariantReadCountLike


class FilterLFN(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.FilterLFN"
    }

    # Input file
    __input_file_sortedinfo = "sortedinfo"
    # Input table
    __input_table_run = "Run"
    __input_table_marker = "Marker"
    __input_table_sample = "Sample"
    __input_table_variant_read_count = "VariantReadCount"
    # Output table
    __output_table_filter_lfn = "FilterLFN"

    def specify_input_file(self):
        return[
            "sortedinfo", "params", "cutoff_specific"
        ]

    def specify_input_table(self):
        return [
            FilterLFN.__input_table_marker,
            FilterLFN.__input_table_run,
            FilterLFN.__input_table_sample,
            FilterLFN.__input_table_variant_read_count,
        ]

    def specify_output_table(self):
        return [
            FilterLFN.__output_table_filter_lfn,
        ]

    def specify_params(self):
        return {
            "lfn_variant_cutoff": "float",
            "lfn_variant_specific_cutoff": "str",
            "lfn_variant_replicate_cutoff": "float",
            "lfn_variant_replicate_specific_cutoff": "str",
            "lfn_sample_replicate_cutoff": "required|float",
            "lfn_read_count_cutoff": "required|float",
        }

    def run(self):

        session = self.session
        engine = session._session().get_bind()

        ############################################################################################

        #
        # Wrapper inputs, outputs and parameters
        #
        ############################################################################################

        # Input file output
        fasta_info_tsv = self.input_file(FilterLFN.__input_file_sortedinfo)

        #
        # Input table models
        input_variant_read_count_model = self.input_table(
            FilterLFN.__input_table_variant_read_count)
        #
        # Output table models
        output_filter_lfn_model = self.output_table(
            FilterLFN.__output_table_filter_lfn)
        #
        # Options
        lfn_variant_cutoff = self.option("lfn_variant_cutoff")
        lfn_variant_specific_cutoff = self.option("lfn_variant_specific_cutoff")
        lfn_variant_replicate_cutoff = self.option("lfn_variant_replicate_cutoff")
        lfn_variant_replicate_specific_cutoff = self.option("lfn_variant_replicate_specific_cutoff")
        lfn_sample_replicate_cutoff = self.option("lfn_sample_replicate_cutoff")
        lfn_read_count_cutoff = self.option("lfn_read_count_cutoff")

        ############################################################################################
        #
        # 1. Read sortedinfo to get run_id, marker_id, sample_id, replicate for current analysis
        # 2. Delete marker_name/run_name/sample/replicate from variant_read_count_model
        # 3. Get nijk_df input
        #
        ############################################################################################

        sample_info_tsv_obj = FileSampleInformation(tsv_path=fasta_info_tsv)

        sample_info_tsv_obj.delete_from_db(
            engine=engine, variant_read_count_like_model=output_filter_lfn_model)

        variant_read_count_df = sample_info_tsv_obj.get_nijk_df(
            variant_read_count_like_model=input_variant_read_count_model, engine=engine, filter_id=None)

        lfn_variant_specific_cutoff_df = None
        if (not (lfn_variant_cutoff is None)) and pathlib.Path(lfn_variant_specific_cutoff).stat().st_size > 0:
            lfn_variant_specific_cutoff_df = FileCutoffSpecific(lfn_variant_specific_cutoff).to_identifier_df(engine=engine, is_lfn_variant_replicate=False)

        lfn_variant_replicate_specific_cutoff_df = None
        if (not (lfn_variant_replicate_cutoff is None)) and pathlib.Path(lfn_variant_replicate_specific_cutoff).stat().st_size > 0:
            lfn_variant_replicate_specific_cutoff_df = FileCutoffSpecific(lfn_variant_replicate_specific_cutoff).to_identifier_df(engine=engine, is_lfn_variant_replicate=True)

        ############################################################################################
        #
        # Create filter object and run_name
        #
        ############################################################################################

        variant_read_count_delete_df = RunnerFilterLFN(variant_read_count_df).get_variant_read_count_delete_df(
            lfn_variant_cutoff=lfn_variant_cutoff,
            lfn_variant_specific_cutoff=lfn_variant_specific_cutoff_df,
            lfn_variant_replicate_cutoff=lfn_variant_replicate_cutoff,
            lfn_variant_replicate_specific_cutoff=lfn_variant_replicate_specific_cutoff_df,
            lfn_sample_replicate_cutoff=lfn_sample_replicate_cutoff,
            lfn_read_count_cutoff=lfn_read_count_cutoff)

        DataframeVariantReadCountLike(variant_read_count_delete_df).to_sql(
            engine=engine, variant_read_count_like_model=output_filter_lfn_model)

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
