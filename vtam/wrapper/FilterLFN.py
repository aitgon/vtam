from vtam.utils.FilterLFNRunner import FilterLFNrunner
from vtam.utils.Logger import Logger
from vtam.utils.SampleInformationFile import SampleInformationFile
from vtam.utils.VTAMexception import VTAMexception
from wopmars.models.ToolWrapper import ToolWrapper

import sys

from vtam.utils.VariantReadCountLikeTable import VariantReadCountLikeTable


class FilterLFN(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.FilterLFN"
    }

    # Input file
    __input_file_readinfo = "readinfo"
    # Input table
    __input_table_run = "Run"
    __input_table_marker = "Marker"
    __input_table_biosample = "Biosample"
    __input_table_variant_read_count = "VariantReadCount"
    # Output table
    __output_table_filter_lfn = "FilterLFN"

    def specify_input_file(self):
        return[
            FilterLFN.__input_file_readinfo,
        ]

    def specify_input_table(self):
        return [
            FilterLFN.__input_table_marker,
            FilterLFN.__input_table_run,
            FilterLFN.__input_table_biosample,
            FilterLFN.__input_table_variant_read_count,
        ]

    def specify_output_table(self):
        return [
            FilterLFN.__output_table_filter_lfn,
        ]

    def specify_params(self):
        return {
            "lfn_variant_threshold": "float",
            "lfn_variant_replicate_threshold": "float",
            "lfn_biosample_replicate_threshold": "required|float",
            "lfn_read_count_threshold": "required|float",
        }

    def run(self):

        session = self.session
        engine = session._session().get_bind()

        ################################################################################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ################################################################################################################

        # Input file output
        fasta_info_tsv = self.input_file(FilterLFN.__input_file_readinfo)
        # Add FilterLFNthresholdspecific
        # input_file_threshold_specific = self.input_file(FilterLFNthresholdspecific.__input_file_threshold_specific)
        #
        # Input table models
        VariantReadCount = self.input_table(FilterLFN.__input_table_variant_read_count)
        #
        # Output table models
        output_filter_lfn_model = self.output_table(FilterLFN.__output_table_filter_lfn)
        #
        # Options
        lfn_variant_threshold = self.option("lfn_variant_threshold")
        lfn_variant_replicate_threshold = self.option("lfn_variant_replicate_threshold")
        lfn_biosample_replicate_threshold = self.option("lfn_biosample_replicate_threshold")
        lfn_read_count_threshold = self.option("lfn_read_count_threshold")

        ################################################################################################################
        #
        # 1. Read readinfo to get run_id, marker_id, biosample_id, replicate for current analysis
        # 2. Delete marker/run/biosample/replicate from variant_read_count_model
        # 3. Select marker/run/biosample/replicate from variant_read_count_model
        # 4. Apply filters
        # 5. Write filters to variant_filter_lfn_model
        #
        ################################################################################################################

        ################################################################################################################
        #
        # 1. Read readinfo to get run_id, marker_id, biosample_id, replicate for current analysis
        #
        ################################################################################################################

        sample_info_tsv_obj = SampleInformationFile(tsv_path=fasta_info_tsv)

        ################################################################################################################
        #
        # 2. Delete marker/run/biosample/replicate from variant_read_count_model
        #
        ################################################################################################################

        variant_read_count_like_utils = VariantReadCountLikeTable( variant_read_count_like_model=output_filter_lfn_model, engine=engine)
        sample_record_list = sample_info_tsv_obj.to_identifier_df(engine=engine).to_dict('records')
        variant_read_count_like_utils.delete_from_db(sample_record_list=sample_record_list)

        ################################################################################################################
        #
        # 3. Select marker/run/biosample/replicate from variant_read_count_model
        #
        ################################################################################################################

        variant_read_count_df = sample_info_tsv_obj.get_variant_read_count_df(
            variant_read_count_like_model=VariantReadCount, engine=engine, filter_id=None)

        ################################################################################################################
        #
        # Create filter object and run
        #
        ################################################################################################################

        lfn_filter_runner = FilterLFNrunner(variant_read_count_df)
        #
        Logger.instance().info("Launching LFN filter:")

        filter_output_df = lfn_filter_runner.run(lfn_variant_threshold, lfn_variant_replicate_threshold,
                                                 lfn_biosample_replicate_threshold, lfn_read_count_threshold)

        ################################################################################################################
        #
        # Write to DB
        #
        ################################################################################################################

        # filter_output_df = lfn_filter_runner.variant_read_count_filter_delete_df

        record_list = VariantReadCountLikeTable.filter_delete_df_to_dict(filter_output_df)

        with engine.connect() as conn:

            # Insert new instances
            conn.execute(output_filter_lfn_model.__table__.insert(), record_list)

        ################################################################################################################
        #
        # Touch output tables, to update modification date
        #
        ################################################################################################################

        for output_table_i in self.specify_output_table():

            declarative_meta_i = self.output_table(output_table_i)
            obj = session.query(declarative_meta_i).order_by(declarative_meta_i.id.desc()).first()
            session.query(declarative_meta_i).filter_by(id=obj.id).update({'id': obj.id})
            session.commit()

        ##########################################################
        #
        # Exit vtam if all variants deleted
        #
        ##########################################################

        if filter_output_df.filter_delete.sum() == filter_output_df.shape[0]:
            Logger.instance().warning(VTAMexception("This filter has deleted all the variants: {}. "
                                                    "The analysis will stop here.".format(self.__class__.__name__)))
            sys.exit(0)


