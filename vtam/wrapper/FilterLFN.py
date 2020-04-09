import sqlalchemy
from sqlalchemy import bindparam

from vtam.utils.FilterLFNrunner import FilterLFNrunner
from vtam.utils.Logger import Logger
from vtam.utils.SampleInformationUtils import FastaInformationTSV
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
        run_model = self.input_table(FilterLFN.__input_table_run)
        marker_model = self.input_table(FilterLFN.__input_table_marker)
        biosample_model = self.input_table(FilterLFN.__input_table_biosample)
        input_variant_read_count_model = self.input_table(FilterLFN.__input_table_variant_read_count)
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

        fasta_info_tsv = FastaInformationTSV(engine=engine, fasta_info_tsv=fasta_info_tsv)

        ################################################################################################################
        #
        # 2. Delete marker/run/biosample/replicate from variant_read_count_model
        #
        ################################################################################################################

        variant_read_count_like_utils = VariantReadCountLikeTable(variant_read_count_like_model=output_filter_lfn_model,
                                                                  engine=engine)
        variant_read_count_like_utils.delete_from_db(sample_record_list=fasta_info_tsv.sample_record_list)

        ################################################################################################################
        #
        # 3. Select marker/run/biosample/replicate from variant_read_count_model
        #
        ################################################################################################################

        variant_read_count_df = fasta_info_tsv.get_variant_read_count_df(
            variant_read_count_like_model=input_variant_read_count_model, filter_id=None)

        ################################################################################################################
        #
        #
        ################################################################################################################

        #
        lfn_filter_runner = FilterLFNrunner(variant_read_count_df)
        #
        Logger.instance().info("Launching LFN filter:")

        ################################################################################################################
        #
        # Filter 2: f2_f4_lfn_delete_variant
        # Or
        # Filter  3: f3_f5_lfn_delete_variant_replicate
        #
        ################################################################################################################

        if lfn_variant_replicate_threshold is None:  # run lfn_variant
            # lfn_filter_runner.f2_f4_lfn_delete_variant(lfn_variant_threshold)
            lfn_filter_runner.mark_delete_lfn_per_Ni_or_Nik_or_Njk(lfn_denominator='N_i', threshold=lfn_variant_threshold)
        else:  # run lfn_variant_replicate
            # lfn_filter_runner.f3_f5_lfn_delete_variant_replicate(lfn_variant_replicate_threshold)
            lfn_filter_runner.mark_delete_lfn_per_Ni_or_Nik_or_Njk(lfn_denominator='N_ik', threshold=lfn_variant_replicate_threshold)

        ################################################################################################################
        #
        # Filter 6:  f6_lfn_delete_biosample_replicate_delete
        #
        ################################################################################################################

        lfn_filter_runner.mark_delete_lfn_per_Ni_or_Nik_or_Njk(lfn_denominator='N_jk', threshold=lfn_biosample_replicate_threshold)

        ################################################################################################################
        #
        # Filter  7:mark_delete_lfn_absolute_read_count
        #
        ################################################################################################################

        lfn_filter_runner.mark_delete_lfn_absolute_read_count(lfn_read_count_threshold)

        ################################################################################################################
        #
        # Filter 8:mark_delete_lfn_do_not_pass_all_filters
        #
        ################################################################################################################

        lfn_filter_runner.mark_delete_lfn_do_not_pass_all_filters()

        ################################################################################################################
        #
        # Write to DB
        #
        ################################################################################################################

        filter_output_df = lfn_filter_runner.variant_read_count_filter_delete_df

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


