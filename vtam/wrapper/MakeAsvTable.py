import inspect
from wopmars.models.ToolWrapper import ToolWrapper
from sqlalchemy import select
import pandas

from vtam.utils.FastaInformation import FastaInformation
from vtam.utils.MakeAsvTableRunner import MakeAsvTableRunner
from vtam.utils.OptionManager import OptionManager
from vtam.utils.Logger import Logger
from vtam.utils.TaxAssignUtilities import f04_1_tax_id_to_taxonomy_lineage, f01_taxonomy_sqlite_to_df
from vtam.utils.constants import rank_hierarchy_asv_table


class MakeAsvTable(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.MakeAsvTable"}

    # Input file
    __input_file_fastainfo = "fastainfo"
    __input_file_taxonomy = "taxonomy"
    # Input table
    __input_table_marker = "Marker"
    __input_table_run = "Run"
    __input_table_biosample = "Biosample"
    __input_table_replicate = "Replicate"
    __input_table_filter_chimera_borderline = "FilterChimeraBorderline"
    __input_table_filter_codon_stop = "FilterCodonStop"
    __input_table_variant = "Variant"
    __input_table_tax_assign = "TaxAssign"

    # Output table
    __output_table_asv = "ASVTable"

    def specify_input_file(self):
        return[
            MakeAsvTable.__input_file_fastainfo,
            MakeAsvTable.__input_file_taxonomy,
        ]

    def specify_input_table(self):
        return [
            MakeAsvTable.__input_table_marker,
            MakeAsvTable.__input_table_run,
            MakeAsvTable.__input_table_biosample,
            MakeAsvTable.__input_table_replicate,
            MakeAsvTable.__input_table_variant,
            MakeAsvTable.__input_table_filter_chimera_borderline,
            MakeAsvTable.__input_table_filter_codon_stop,
            MakeAsvTable.__input_table_tax_assign,
        ]

    def specify_output_file(self):
        return[
            MakeAsvTable.__output_table_asv,

        ]

    def specify_params(self):
        return {
            "foo": "int",
            # "log_verbosity": "int",
            # "log_file": "str"
        }

    def run(self):
        session = self.session
        engine = session._session().get_bind()

        # if not self.option("log_verbosity") is None:
        #     OptionManager.instance()['log_verbosity'] = int(self.option("log_verbosity"))
        #     OptionManager.instance()['log_file'] = str(self.option("log_file"))

        ##########################################################
        #
        # 1. Wrapper inputs, outputs and parameters
        #
        ##########################################################
        #
        # Input file output
        input_file_fastainfo = self.input_file(MakeAsvTable.__input_file_fastainfo)
        input_file_taxonomy = self.input_file(MakeAsvTable.__input_file_taxonomy)
        #
        # Input table models
        marker_model = self.input_table(MakeAsvTable.__input_table_marker)
        run_model = self.input_table(MakeAsvTable.__input_table_run)
        biosample_model = self.input_table(MakeAsvTable.__input_table_biosample)
        replicate_model = self.input_table(MakeAsvTable.__input_table_replicate)
        filter_chimera_borderline_model = self.input_table(MakeAsvTable.__input_table_filter_chimera_borderline)
        filter_codon_stop_model = self.input_table(MakeAsvTable.__input_table_filter_codon_stop)
        variant_model = self.input_table(MakeAsvTable.__input_table_variant)
        tax_assign_model = self.input_table(MakeAsvTable.__input_table_tax_assign)
        # Output table models
        asv_table_tsv_path = self.output_file(MakeAsvTable.__output_table_asv)
        #
        # Options

        ##########################################################
        #
        # 1. Read fastainfo to get run_id, marker_id, biosample_id, replicate_id for current analysis
        # Compute variant_read_count_df
        #
        ##########################################################

        asv_table_runner = MakeAsvTableRunner(engine, input_file_fastainfo, run_model, marker_model, biosample_model, replicate_model, filter_chimera_borderline_model,
                                              filter_codon_stop_model, variant_model, tax_assign_model, input_file_taxonomy)
        asv_df_final = asv_table_runner.run()

        asv_df_final.to_csv(asv_table_tsv_path, sep='\t', index=False, header=True)
