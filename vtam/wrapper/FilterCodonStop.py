from vtam.utils.RunnerFilterCodonStop import RunnerFilterCodonStop
from vtam.utils.Logger import Logger
from vtam.utils.FileSampleInformation import FileSampleInformation
from vtam.utils.VTAMexception import VTAMexception
from vtam.utils.DataframeVariantReadCountLike import DataframeVariantReadCountLike
from wopmars.models.ToolWrapper import ToolWrapper
import sys


class FilterCodonStop(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.FilterCodonStop"
    }

    # Input file
    __input_file_sortedinfo = "sortedinfo"
    # Input table
    __input_table_marker = "Marker"
    __input_table_run = "Run"
    __input_table_sample = "Sample"
    __input_table_filter_indel = "FilterIndel"
    __input_table_Variant = "Variant"
    # Output table
    __output_table_filter_codon_stop = "FilterCodonStop"

    def specify_input_file(self):
        return[
            "sortedinfo", "params",
        ]

    def specify_input_table(self):
        return [
            FilterCodonStop.__input_table_marker,
            FilterCodonStop.__input_table_run,
            FilterCodonStop.__input_table_sample,
            FilterCodonStop.__input_table_filter_indel,
            FilterCodonStop.__input_table_Variant,
        ]

    def specify_output_table(self):
        return [
            FilterCodonStop.__output_table_filter_codon_stop,
        ]

    def specify_params(self):
        return {
            "genetic_code": "int",
            "skip_filter_codon_stop": "int",
        }

    def run(self):
        session = self.session
        engine = session._session().get_bind()

        ##########################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ##########################################################
        #
        # Input file output
        fasta_info_tsv = self.input_file(FilterCodonStop.__input_file_sortedinfo)
        #
        # Input table models
        input_filter_indel_model = self.input_table(
            FilterCodonStop.__input_table_filter_indel)
        #
        # Options
        genetic_code = int(self.option("genetic_code"))
        skip_filter_codon_stop = bool(int(self.option("skip_filter_codon_stop")))
        #
        # Output table models
        output_filter_codon_stop_model = self.output_table(
            FilterCodonStop.__output_table_filter_codon_stop)

        #######################################################################
        #
        # 1. Read sortedinfo to get run_id, marker_id, sample_id, replicate for current analysis
        # 2. Delete marker_name/run_name/sample/replicate from variant_read_count_model
        # 3. Get nijk_df input
        #
        #######################################################################

        sample_info_tsv_obj = FileSampleInformation(tsv_path=fasta_info_tsv)

        sample_info_tsv_obj.delete_from_db(
            engine=engine, variant_read_count_like_model=output_filter_codon_stop_model)

        variant_read_count_df = sample_info_tsv_obj.get_nijk_df(
            variant_read_count_like_model=input_filter_indel_model, engine=engine, filter_id=None)

        #######################################################################
        #
        # 4. Run Filter
        #
        #######################################################################

        variant_df = sample_info_tsv_obj.get_variant_df(
            variant_read_count_like_model=input_filter_indel_model, engine=engine)
        variant_read_count_delete_df = RunnerFilterCodonStop(
            variant_read_count_df=variant_read_count_df).get_variant_read_count_delete_df(
            variant_df=variant_df,
            genetic_code=genetic_code,
            skip_filter_codon_stop=skip_filter_codon_stop)

        #######################################################################
        #
        # 5. Write to DB
        # 6. Touch output tables, to update modification date
        # 7. Exit vtam if all variants delete
        #
        #######################################################################

        DataframeVariantReadCountLike(variant_read_count_delete_df).to_sql(
            engine=engine, variant_read_count_like_model=output_filter_codon_stop_model)

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
