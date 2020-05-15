from vtam.utils.FilterCodonStopRunner2 import FilterCodonStopRunner2
from vtam.utils.SampleInformationFile import SampleInformationFile
from vtam.utils.VariantReadCountLikeTable import VariantReadCountLikeTable
from vtam.utils.Logger import Logger
from vtam.utils.VTAMexception import VTAMexception
from wopmars.models.ToolWrapper import ToolWrapper

import sys


class FilterCodonStop(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.FilterCodonStop"
    }

    # Input file
    __input_file_readinfo = "readinfo"
    # Input table
    __input_table_marker = "Marker"
    __input_table_run = "Run"
    __input_table_biosample = "Biosample"
    __input_table_filter_indel = "FilterIndel"
    __input_table_Variant = "Variant"
    # Output table
    __output_table_filter_codon_stop = "FilterCodonStop"

    def specify_input_file(self):
        return[
            FilterCodonStop.__input_file_readinfo,
        ]

    def specify_input_table(self):
        return [
            FilterCodonStop.__input_table_marker,
            FilterCodonStop.__input_table_run,
            FilterCodonStop.__input_table_biosample,
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
        fasta_info_tsv = self.input_file(FilterCodonStop.__input_file_readinfo)
        #
        # Input table models
        Variant = self.input_table(FilterCodonStop.__input_table_Variant)
        input_filter_indel_model = self.input_table(FilterCodonStop.__input_table_filter_indel)
        #
        # Options
        genetic_code = int(self.option("genetic_code"))
        skip_filter_codon_stop = bool(int(self.option("genetic_code")))
        #
        # Output table models
        output_filter_codon_stop_model = self.output_table(FilterCodonStop.__output_table_filter_codon_stop)

        ################################################################################################################
        #
        # 1. Read readinfo to get run_id, marker_id, biosample_id, replicate for current analysis
        #
        ################################################################################################################

        # fasta_info_tsv = FastaInformationTSV(engine=engine, fasta_info_tsv=fasta_info_tsv)
        sample_info_tsv_obj = SampleInformationFile(tsv_path=fasta_info_tsv)

        ################################################################################################################
        #
        # 2. Delete marker/run/biosample/replicate from variant_read_count_model
        #
        ################################################################################################################

        variant_read_count_like_utils = VariantReadCountLikeTable(variant_read_count_like_model=output_filter_codon_stop_model, engine=engine)
        sample_record_list = sample_info_tsv_obj.to_identifier_df(engine=engine).to_dict('records')
        variant_read_count_like_utils.delete_from_db(sample_record_list=sample_record_list)

        ################################################################################################################
        #
        # variant_read_count_input_df
        #
        ################################################################################################################

        variant_read_count_df = sample_info_tsv_obj.get_variant_read_count_df(
            variant_read_count_like_model=input_filter_indel_model, engine=engine, filter_id=None)
        variant_df = sample_info_tsv_obj.get_variant_df(variant_read_count_like_model=input_filter_indel_model, engine=engine)

        ##########################################################
        #
        # 4. Run Filter or not according to skip_filter_codon_stop
        #
        ##########################################################

        if skip_filter_codon_stop:  # do not run filter

            filter_output_df = variant_read_count_df.copy()
            filter_output_df['filter_delete'] = False

        else:  # run filter
            filter_codon_stop_runner_obj = FilterCodonStopRunner2(code=genetic_code)
            variant_has_stop_codon_df = filter_codon_stop_runner_obj.annotate_stop_codon_count(variant_df)
            variants_with_stop_codons_list = variant_has_stop_codon_df.id[
                variant_has_stop_codon_df['has_stop_codon'] == 1].tolist()

            filter_output_df = variant_read_count_df.copy()
            filter_output_df['filter_delete'] = False
            filter_output_df.loc[filter_output_df.variant_id.isin(filter_output_df), 'filter_delete'] = True
            # filter_output_df = f14_filter_codon_stop(variant_read_count_df, variant_df, genetic_code)

        ##########################################################
        #
        # Write to DB
        #
        ##########################################################

        record_list = VariantReadCountLikeTable.filter_delete_df_to_dict(filter_output_df)
        with engine.connect() as conn:

            # Insert new instances
            conn.execute(output_filter_codon_stop_model.__table__.insert(), record_list)

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
