from wopmars.models.ToolWrapper import ToolWrapper

from vtam.utils.RunnerAsvTable import RunnerAsvTable
from vtam.utils.FileKnownOccurrences import FileKnownOccurrences
from vtam.utils.FileSampleInformation import FileSampleInformation
from vtam.models.FilterCodonStop import FilterCodonStop


class MakeAsvTable(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.MakeAsvTable"}

    # Input file
    __input_file_sortedinfo = "sortedinfo"
    __input_file_known_occurrences = "known_occurrences"
    # Input table
    __input_table_marker = "Marker"
    __input_table_run = "Run"
    __input_table_sample = "Sample"
    __input_table_filter_chimera_borderline = "FilterChimeraBorderline"
    __input_table_filter_codon_stop = "FilterCodonStop"
    __input_table_variant = "Variant"
    __input_table_tax_assign = "TaxAssign"

    # Output table
    __output_table_asv = "ASVTable"

    def specify_input_file(self):
        return[
            MakeAsvTable.__input_file_sortedinfo,
        ]

    def specify_input_table(self):
        return [
            MakeAsvTable.__input_table_marker,
            MakeAsvTable.__input_table_run,
            MakeAsvTable.__input_table_sample,
            MakeAsvTable.__input_table_variant,
            MakeAsvTable.__input_table_filter_chimera_borderline,
            MakeAsvTable.__input_table_filter_codon_stop,
        ]

    def specify_output_file(self):
        return[
            MakeAsvTable.__output_table_asv,

        ]

    def specify_params(self):
        return {
            "foo": "int",
            "cluster_identity": "float",
            "known_occurrences": "str",
        }

    def run(self):
        session = self.session
        engine = session._session().get_bind()

        #######################################################################
        #
        # 1. Wrapper inputs, outputs and parameters
        #
        #######################################################################

        # Input file
        fasta_info_tsv = self.input_file(MakeAsvTable.__input_file_sortedinfo)

        # Output file
        asvtable_tsv_path = self.output_file(MakeAsvTable.__output_table_asv)
        #
        # Options
        cluster_identity = float(self.option("cluster_identity"))
        known_occurrences_tsv = str(self.option("known_occurrences"))

        #######################################################################
        #
        # Read sortedinfo to get run_id, marker_id, sample_id, replicate for current analysis
        # Compute variant_read_count_input_df and other dfs for the asv_table_runner
        #
        #######################################################################

        sample_info_tsv_obj = FileSampleInformation(tsv_path=fasta_info_tsv)

        variant_read_count_df = sample_info_tsv_obj.get_nijk_df(
            FilterCodonStop, engine=engine)

        ############################################################################################
        #
        # FileKnownOccurrences
        #
        ############################################################################################

        if known_occurrences_tsv == 'None' or known_occurrences_tsv is None:
            known_occurrences_df = None
        else:
            known_occurrences_df = FileKnownOccurrences(known_occurrences_tsv).to_identifier_df(engine)
            known_occurrences_df = known_occurrences_df.loc[
                (known_occurrences_df.mock == 1) & (known_occurrences_df.action == 'keep'), ]

        #######################################################################
        #
        # Compute variant_to_chimera_borderline_df
        #
        #######################################################################

        sample_list = sample_info_tsv_obj.read_tsv_into_df()['sample'].drop_duplicates(keep='first').tolist()
        asvtable_runner = RunnerAsvTable(variant_read_count_df=variant_read_count_df,
                                         engine=engine, sample_list=sample_list,
                                         cluster_identity=cluster_identity,
                                         known_occurrences_df=known_occurrences_df)
        asvtable_runner.to_tsv(asvtable_tsv_path)
