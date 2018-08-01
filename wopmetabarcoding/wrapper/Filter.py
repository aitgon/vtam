import errno
from wopmars.utils.Logger import Logger
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger

from wopmetabarcoding.utils.constants import tempdir
from wopmetabarcoding.wrapper.FilterUtilities import Variant2Sample2Replicate2Count
from sqlalchemy import select
import pandas, os, pickle


class Filter(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.Filter"
    }

    # Input tables:
    __input_table_marker = "Marker"
    __input_table_variant = "Variant"
    __input_table_variant_read_count = "VariantReadCount"
    __input_table_replicate = "Replicate"
    __output_table_passed_variant = "PassedVariant"
    # Input file
    __input_sortreads_samplecount_csv = "SortReads_sample_counts"
    __input_cutoff_file = "file_cutoff"
    __input_genetic_code_file = "genetic_code"
    # Output file
    __output_marker_variant_path = "marker_variant_path"


    def specify_input_table(self):
        return [
            Filter.__input_table_marker,
            Filter.__input_table_variant,
            Filter.__input_table_variant_read_count,
            Filter.__input_table_replicate
        ]

    def specify_input_file(self):
        return[
            # Filter.__input_sortreads_samplecount_csv,
            Filter.__input_cutoff_file,
            Filter.__input_genetic_code_file
        ]


    def specify_output_table(self):
        return [
            Filter.__output_table_passed_variant,
        ]

    def specify_output_file(self):
        return[
            Filter.__output_marker_variant_path
        ]

    def specify_params(self):
        return {
            "filter_output_dir": "str",
            "lfn_per_replicate_threshold": "float",
            "lfn_per_variant_threshold": "float",
            "lfn_per_replicate_series_threshold": "float",
            "lfn_read_count_threshold": "float",
            "repeatability_threshold": "int",
            "pcr_error_var_prop": "float",
            "renkonen_number_of_replicate": "int",
            "renkonen_renkonen_tail": "float"
        }

    # Todo: Needs output fields, eg table with paths to marker data frame pickles

    @staticmethod
    def codon_stop_dataframe(genetic_code_csv):
        # df_codon_stop_columns = ['codon']
        # df_codon_stop = pandas.DataFrame(columns=df_codon_stop_columns)
        df_genetic_code_columns = ['genetic_code', 'codon']
        df_genetic_code = pandas.DataFrame(columns=df_genetic_code_columns)
        i = 1
        with open(genetic_code_csv, 'r') as fin:
            lines = []
            for line in fin:
                if i <= 3:
                    lines.append(line)
                if i >= 4:
                    base1 = lines[0]
                    base2 = lines[1]
                    base3 = lines[2]
                    for linebis in fin:
                        linebis = linebis.split(';')
                        j = 0
                        sequence = linebis[0]
                        for j in range(len(sequence)):
                            if sequence[j] == '*':
                                codon = base1[j] + base2[j] + base3[j]
                                df_genetic_code.loc[len(df_genetic_code)] = [linebis[1], codon]
                                j += 1
                i += 1
        return df_genetic_code

            # base1 = lines[0]
            # base2 = lines[1]
            # base3 = lines[2]
            # j = 0
            # for j in range(len(base1)):
            #     if base1[j] == ";":
            #         break
            #     codon = base1[j] + base2[j] + base3[j]
            #     df_codon_stop.loc[len(df_codon_stop)] = [codon]
            #     j += 1

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        conn = engine.connect()
        #
        # # Input file path
        # # sortreads_samplecount = self.input_file(Filter.__input_sortreads_samplecount_csv)
        cutoff_file_tsv = self.input_file(Filter.__input_cutoff_file)
        genetic_code_tsv = self.input_file(Filter.__input_genetic_code_file)
        #
        # Input table models
        marker_model = self.input_table(Filter.__input_table_marker)
        variant_model = self.input_table(Filter.__input_table_variant)
        variant_read_count_model = self.input_table(Filter.__input_table_variant_read_count)
        replicate_model = self.input_table(Filter.__input_table_replicate)
        #
        # Output table models
        dropped_variant_model = self.output_table(Filter.__output_table_passed_variant)
        #
        variant2biosample2replicate2count_list = []
        variant_read_count_model_table = variant_read_count_model.__table__
        stmt_marker_id = select([variant_read_count_model_table.c.marker_id]).distinct()
        for row in conn.execute(stmt_marker_id).fetchall():
            marker_id = row[0]
            stmt_variant = select([variant_read_count_model_table.c.variant_id, variant_read_count_model_table.c.biosample_id, variant_read_count_model_table.c.replicate_id, variant_read_count_model_table.c.read_count]).where(
                variant_read_count_model_table.c.marker_id == marker_id)
            for row2 in conn.execute(stmt_variant).fetchall():
                variant2biosample2replicate2count_list.append(row2)
            variant2biosample2replicate2count_df = pandas.DataFrame.from_records(variant2biosample2replicate2count_list, columns=['variant_id', 'biosample_id', 'replicate_id', 'read_count'])
            #
            filter_obj = Variant2Sample2Replicate2Count(variant2biosample2replicate2count_df)
            #
            # Filter parameters
            lfn_per_replicate_threshold = self.option("lfn_per_replicate_threshold")
            lfn_per_variant_threshold = self.option("lfn_per_variant_threshold")
            lfn_per_replicate_series_threshold = self.option("lfn_per_replicate_series_threshold")
            lfn_read_count_threshold = self.option("lfn_read_count_threshold")
            #
            Logger.instance().info("Launching LFN filter:")
            #
            ############################################
            # Filter 1: f1_lfn1_per_replicate
            ############################################
            filter_obj.f1_lfn1_per_replicate(lfn_per_replicate_threshold)
            #
            ############################################
            # Filter 2: f2_lfn2_per_variant
            # Filter 3: f3_lfn2_per_replicate_series
            ############################################
            f2_lfn2_per_variant = False
            if f2_lfn2_per_variant:
                filter_obj.f2_lfn2_per_variant(lfn_per_variant_threshold)
            else:
                filter_obj.f3_lfn2_per_replicate_series(lfn_per_replicate_series_threshold)
            #
            ############################################
            # Filter 4: lfn3_read_count
            ############################################
            lfn3_to_drop_variant_id_list = filter_obj.lfn3_read_count(lfn_read_count_threshold)
            #
            ############################################
            # Filter 5: lfn4_per_variant_with_cutoff
            # Filter 6: lfn4_per_replicate_series_with_cutoff
            ############################################
            if not lfn_per_replicate_series:
                filter_obj.lfn4_per_variant_with_cutoff(cutoff_file_tsv)
            else:
                filter_obj.lfn4_per_replicate_series_with_cutoff(cutoff_file_tsv)
            #
            ############################################
            # Filter 7: Repeatability: min_repln
            # Filter 8: Repeatability: min_replp
            ############################################
            replp = False
            if replp is False:
                filter_obj.min_repln(2)
            else:
                filter_obj.min_replp(3)
            #
            ############################################
            # Filter 9: PCR error
            ############################################
            # Logger.instance().info("Launching PCR error filter:")
            # filter_obj.store_variant_ids_identified_as_pcr_error(biosample_list, sample_replicate_list, marker_id, self.option("pcr_error_var_prop"), False)
            # filter_obj.execute_filters()
            #
            ############################################
            # Filter 10: Chimera
            ############################################
            # Logger.instance().info("Launching store_index_identified_as_chimerenkonen_renkonen_tailra filter:")
            # # TODO Remove this select and get information from df: variant2sample2replicate2count_df
            # select_replicate = select([replicate_model.biosample_name, replicate_model.name]) \
            #     .where(replicate_model.marker_id == marker_id)
            # replicate_obj = engine.execute(select_replicate)
            # replicate_obj_list = replicate_obj.fetchall()
            # filter_obj.store_variant_ids_identified_as_chimera(replicate_obj_list, marker_id, 'sample_replicate')
            # filter_obj.execute_filters()
            #
            ############################################
            # Filter 11: Renkonen
            ############################################
            # Logger.instance().info("Launching store_variant_ids_that_pass_renkonen filter:")
            # filter_obj.store_variant_ids_that_pass_renkonen(self.option("renkonen_number_of_replicate"), self.option("renkonen_renkonen_tail"))
            # filter_obj.execute_filters()
            #
            ############################################
            # Filter 12: Indel
            ############################################
            #     #
            #     # Filter: indel
            #     Logger.instance().info("Launching pseudogene detection with indel filter:")
            #     # TODO check if here it is empty
            #     filter_obj.indel(delete_var=False)
            #     df_codon_stop_per_genetic_code = self.codon_stop_dataframe(genetic_code_tsv)
            #     Logger.instance().info("Launching pseudogene detection with codon stop filter:")
            #     filter_obj.codon_stop(df_codon_stop_per_genetic_code, 2, False)
            #     filter_obj.consensus()
            #     variant2sample2replicate2count_df = filter_obj.filtered_variants()
            #     filtered_variants_list = list(set(variant2sample2replicate2count_df.variant_seq.tolist()))
            #     output_fasta = os.path.join(filter_output_dir, (marker_name + "_variant.fasta"))
            #     dataframe_tsv_path = os.path.join(filter_output_dir, (marker_name + "_variant_info.tsv"))
            #     filter_obj.filter_fasta(filtered_variants_list, output_fasta, False)
            #     variant2sample2replicate2count_df.to_csv(dataframe_tsv_path, sep='\t', encoding='utf-8', index=False)
            #     fout.write(marker_name + "\t" + dataframe_tsv_path + "\t" + output_fasta + "\n")
            # fout.close()
            # # TODO: Log tables for pcr and store_variant_ids_identified_as_chimera filter
