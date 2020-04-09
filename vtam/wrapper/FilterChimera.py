import os
import pathlib
import sys

from sqlalchemy import bindparam

from vtam.utils.FilterChimeraRunner import FilterChimeraRunner
from vtam.utils.SampleInformationUtils import FastaInformationTSV
from vtam.utils.VariantReadCountLikeTable import VariantReadCountLikeTable
from vtam.utils.Logger import Logger
from vtam.utils.PathManager import PathManager
from vtam.utils.VTAMexception import VTAMexception
from wopmars.models.ToolWrapper import ToolWrapper


class FilterChimera(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.FilterChimera"
    }

    # Input file
    __input_file_readinfo = "readinfo"
    # Input table
    __input_table_marker = "Marker"
    __input_table_run = "Run"
    __input_table_biosample = "Biosample"
    __input_table_filter_pcr_error = "FilterPCRerror"
    __input_table_Variant = "Variant"
    # Output table
    __output_table_filter_chimera = "FilterChimera"
    __output_table_filter_chimera_borderline = "FilterChimeraBorderline"


    def specify_input_file(self):
        return[
            FilterChimera.__input_file_readinfo,

        ]

    def specify_input_table(self):
        return [
            FilterChimera.__input_table_marker,
            FilterChimera.__input_table_run,
            FilterChimera.__input_table_biosample,
            FilterChimera.__input_table_filter_pcr_error,
            FilterChimera.__input_table_Variant,
        ]


    def specify_output_table(self):
        return [
            FilterChimera.__output_table_filter_chimera,
            FilterChimera.__output_table_filter_chimera_borderline,
        ]

    def specify_params(self):
        return{
        }

    def run(self):
        session = self.session
        engine = session._session().get_bind()

        this_temp_dir = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__))
        pathlib.Path(this_temp_dir).mkdir(exist_ok=True)

        ################################################################################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ################################################################################################################
        #
        # Input file output
        fasta_info_tsv = self.input_file(FilterChimera.__input_file_readinfo)
        #
        # Input table models
        marker_model = self.input_table(FilterChimera.__input_table_marker)
        run_model = self.input_table(FilterChimera.__input_table_run)
        biosample_model = self.input_table(FilterChimera.__input_table_biosample)
        variant_model = self.input_table(FilterChimera.__input_table_Variant)
        input_filter_pcr_error_model = self.input_table(FilterChimera.__input_table_filter_pcr_error)
        #
        # Output table models
        output_filter_chimera_model = self.output_table(FilterChimera.__output_table_filter_chimera)
        filter_chimera_borderline_model = self.output_table(FilterChimera.__output_table_filter_chimera_borderline)

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

        variant_read_count_like_table_obj = VariantReadCountLikeTable(variant_read_count_like_model=output_filter_chimera_model, engine=engine)
        variant_read_count_like_table_obj.delete_from_db(sample_record_list=fasta_info_tsv.sample_record_list)

        variant_read_count_like_table_borderline_obj = VariantReadCountLikeTable(variant_read_count_like_model=filter_chimera_borderline_model, engine=engine)
        variant_read_count_like_table_borderline_obj.delete_from_db(sample_record_list=fasta_info_tsv.sample_record_list)

        ################################################################################################################        #
        #
        # 3. Select marker/run/biosample/replicate from variant_read_count_model
        #
        ################################################################################################################

        variant_read_count_df = fasta_info_tsv.get_variant_read_count_df(
            variant_read_count_like_model=input_filter_pcr_error_model, filter_id=None)
        variant_df = fasta_info_tsv.get_variant_df(variant_read_count_like_model=input_filter_pcr_error_model,
                                               variant_model=variant_model)

        ################################################################################################################
        #
        # 4. Run Filter
        #
        ################################################################################################################

        filter_chimera_runner = FilterChimeraRunner(variant_df=variant_df, variant_read_count_df=variant_read_count_df)
        filter_output_chimera_df, filter_borderline_output_df = filter_chimera_runner.run(tmp_dir=this_temp_dir)

        ################################################################################################################
        #
        # Write to DB
        #
        ################################################################################################################

        records_chimera_list = VariantReadCountLikeTable.filter_delete_df_to_dict(filter_output_chimera_df)

        with engine.connect() as conn:

            # Insert new instances
            conn.execute(output_filter_chimera_model.__table__.insert(), records_chimera_list)

        records_chimera_borderline_list = VariantReadCountLikeTable.filter_delete_df_to_dict(filter_borderline_output_df)

        with engine.connect() as conn:

            # Insert new instances
            conn.execute(filter_chimera_borderline_model.__table__.insert(), records_chimera_borderline_list)

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

        ################################################################################################################
        #
        # Exit vtam if all variants delete
        #
        ################################################################################################################

        if filter_output_chimera_df.filter_delete.sum() == filter_output_chimera_df.shape[0]:
            Logger.instance().warning(VTAMexception("This filter has deleted all the variants: {}. "
                                                    "The analysis will stop here.".format(self.__class__.__name__)))
            sys.exit(0)
