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
    __input_table_replicate = "Replicate"
    # Input file
    __input_cutoff_file = "file_cutoff"
    __input_genetic_code_file = "genetic_code"
    # Output file
    __output_filtered_dataframe = "filtered_dataframe_path"

    def specify_input_table(self):
        return [
            Filter.__input_table_marker,
            Filter.__input_table_variant,
            Filter.__input_table_replicate
        ]

    def specify_input_file(self):
        return[
            Filter.__input_cutoff_file,
            Filter.__input_genetic_code_file
        ]

    def specify_output_file(self):
        return[
            Filter.__output_filtered_dataframe
        ]

    def specify_params(self):
        return {
            "lfn_per_replicate_threshold": "float",
            "lfn_per_variant_threshold": "float",
            "lfn_per_replicate_series_threshold": "float",
            "lfn_read_count_threshold": "float",
            "repeatability_threshold": "int",
            "pcr_error_var_prop": "float",
            "renkonen_number_of_replicate": "int",
            "renkonen_renkonen_tail": "float"
        }

    # Todo: Needs output fields, eg table with paths to marker data frame pickles

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
        #
        # Input file path
        cutoff_file_tsv = self.input_file(Filter.__input_cutoff_file)
        genetic_code_tsv = self.input_file(Filter.__input_genetic_code_file)
        #
        # Input table models
        marker_model = self.input_table(Filter.__input_table_marker)
        variant_model = self.input_table(Filter.__input_table_variant)
        replicate_model = self.input_table(Filter.__input_table_replicate)
        # Output file path
        filtered_dataframe_path = self.output_file(Filter.__output_filtered_dataframe)
        #
        fout = open(filtered_dataframe_path, 'w')
        marker_select = select([marker_model.id, marker_model.name])
        marker_obj = engine.execute(marker_select)
        for marker in marker_obj:
            file_name = os.path.join(tempdir, marker.name + '_sample_count.tsv')
            # cols "sequence", "replicate", "sample", "sample_replicate", "count"]
            variant2sample2replicate2count_df = pandas.read_csv(file_name, sep='\t')
            biosample_list = list(set(variant2sample2replicate2count_df['sample'].tolist()))
            sample_replicate_list = list(set(variant2sample2replicate2count_df['sample_replicate'].tolist()))
            filter_obj = Variant2Sample2Replicate2Count(variant2sample2replicate2count_df)
            #
            # Filter parameters
            lfn_per_replicate_threshold = self.option("lfn_per_replicate_threshold")
            lfn_per_variant_threshold = self.option("lfn_per_variant_threshold")
            lfn_per_replicate_series_threshold = self.option("lfn_per_replicate_series_threshold")
            lfn_read_count_threshold = self.option("lfn_read_count_threshold")
            #
            Logger.instance().info("Launching LFN filter:")
            #
            # Filter 1: lfn1 (low frequency noise)
            filter_obj.store_index_below_lfn1_per_replicate(lfn_per_replicate_threshold)
            #
            # Filter 2: lfn2
            lfn_per_replicate_series = False
            if not lfn_per_replicate_series:
                filter_obj.store_index_below_lfn2_per_variant(lfn_per_variant_threshold)
            else:
                filter_obj.store_index_below_lfn2_per_replicate_series(lfn_per_replicate_series_threshold)
            #
            # Filter 3: lfn3
            filter_obj.store_index_below_lfn3_read_count(lfn_read_count_threshold)
            #
            # Filter 4: lfn4
            if not lfn_per_replicate_series:
                filter_obj.store_index_below_lfn4_per_variant_with_cutoff(cutoff_file_tsv)
            else:
                filter_obj.store_index_below_lfn4_per_replicate_series_with_cutoff(cutoff_file_tsv)
            #
            # Drop indices, once for all lfn filters
            filter_obj.drop_indices()
            #
            # Filter : replp and repln (Repeatability filters)
            Logger.instance().info("Launching Repeatability variants filters:")
            replp = False
            if replp is False:
                # TODO: Fix this inconsistency. These filters remove directly indexing, while LFN filters just add indices to an internal list
                filter_obj.store_index_below_min_repln(biosample_list, 2)
            else:
                filter_obj.store_index_below_min_replp(biosample_list, 3)
            filter_obj.drop_indices()
            #
            # Filter: PCR
            Logger.instance().info("Launching PCR error filter:")
            filter_obj.store_index_identified_as_pcr_error(biosample_list, sample_replicate_list, marker.id, self.option("pcr_error_var_prop"), False)
            filter_obj.drop_indices()
            #
            # Filter: Chimera
            Logger.instance().info("Launching store_index_identified_as_chimerenkonen_renkonen_tailra filter:")
            # TODO Remove this select and get information from df: variant2sample2replicate2count_df
            select_replicate = select([replicate_model.biosample_name, replicate_model.name]) \
                .where(replicate_model.marker_id == marker.id)
            replicate_obj = engine.execute(select_replicate)
            replicate_obj_list = replicate_obj.fetchall()
            filter_obj.store_index_identified_as_chimera(replicate_obj_list, marker.id, 'sample_replicate')
            filter_obj.drop_indices()
            #
            # Filter: replp
            if replp is False:
                filter_obj.store_index_below_min_repln(biosample_list, 1)
            else:
                filter_obj.store_index_below_min_replp(biosample_list, 3)
            filter_obj.drop_indices()
            #
            # Filter: Renkonen
            Logger.instance().info("Launching store_index_that_fails_renkonen filter:")
            filter_obj.store_index_that_fails_renkonen(self.option("renkonen_number_of_replicate"), self.option("renkonen_renkonen_tail"))
            filter_obj.drop_indices()
            #
            # Filter: replp
            if replp is False:
                filter_obj.store_index_below_min_repln(biosample_list, 1)
            else:
                filter_obj.store_index_below_min_replp(biosample_list, 3)
            filter_obj.drop_indices()
            #
            # Filter: indel
            Logger.instance().info("Launching pseudogene detection with indel filter:")
            filter_obj.indel(False)
            df_codon_stop_per_genetic_code = self.codon_stop_dataframe(genetic_code_tsv)
            Logger.instance().info("Launching pseudogene detection with codon stop filter:")
            filter_obj.codon_stop(df_codon_stop_per_genetic_code, 2, False)
            filter_obj.consensus()
            variant2sample2replicate2count_df = filter_obj.filtered_variants()
            filtered_variants_list = list(set(variant2sample2replicate2count_df.sequence.tolist()))
            output_fasta = os.path.join(tempdir, (marker.name + "_filtered_variants.fasta"))
            dataframe_pickle_path = os.path.join(tempdir, (marker.name + "_filtered_dataframe.pkl"))
            filter_obj.filter_fasta(filtered_variants_list, output_fasta, False)
            pickle.dump(variant2sample2replicate2count_df, open(dataframe_pickle_path, "wb"))
            fout.write(marker.name + "\t" + dataframe_pickle_path + "\t" + output_fasta + "\n")
        fout.close()
            # TODO: Log tables for pcr and store_index_identified_as_chimera filter
