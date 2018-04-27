from wopmars.utils.Logger import Logger
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger

from wopmetabarcoding.utils.constants import tempdir
from wopmetabarcoding.wrapper.FilterUtilities import Variant2Sample2Replicate2Count
from sqlalchemy import select
import pandas, os


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
        #
        marker_select = select([marker_model.id, marker_model.name])
        marker_obj = engine.execute(marker_select)
        for marker in marker_obj:
            file_name = os.path.join(tempdir, marker.name + '_sample_count.tsv')
            variant2sample2replicate2count_df = pandas.read_csv(file_name, sep='\t')
            filter = Variant2Sample2Replicate2Count(variant2sample2replicate2count_df)
            lfn_per_replicate_threshold = 0.001
            lfn_per_variant_threshold = 0.001
            lfn_per_replicate_series_threshold = 0.001
            lfn_read_count_threshold = 1
            #
            Logger.instance().info("Launching LFN filter:")
            filter.lfn1_per_replicate(lfn_per_replicate_threshold)
            lfn_per_replicate_series = False
            if not lfn_per_replicate_series:
                filter.lfn2_per_variant(lfn_per_variant_threshold)
            else:
                filter.lfn2_per_replicate_series(lfn_per_replicate_series_threshold)
            filter.lfn3_read_count(lfn_read_count_threshold)
            if not lfn_per_replicate_series:
                filter.lfn4_per_variant_with_cutoff(cutoff_file_tsv)
            else:
                filter.lfn4_per_replicate_series_with_cutoff(cutoff_file_tsv)
            #
            filter.drop_indices()
            # min_repln(variant2sample2replicate2count_df, 2)
            Logger.instance().info("Launching Repeatability variants filters:")
            replp = False
            if replp is False:
                filter.min_repln(engine, replicate_model, marker.id, 2)
            else:
                filter.min_replp(engine, replicate_model, marker.id, 3)
            Logger.instance().info("Launching PCR error filter:")
            filter.pcr_error(engine, replicate_model, variant_model, marker.id, 0.5, False)
            Logger.instance().info("Launching chimera filter:")
            select_replicate = select([replicate_model.biosample_name, replicate_model.name]) \
                .where(replicate_model.marker_id == marker.id)
            replicate_obj = engine.execute(select_replicate)
            replicate_obj_list = replicate_obj.fetchall()
            filter.chimera(replicate_obj_list, marker.id, 'sample_replicate')
            if replp is False:
                filter.min_repln(engine, replicate_model, marker.id, 1)
            else:
                filter.min_replp(engine, replicate_model, marker.id, 3)
            filter.renkonen(3, 0.1)
            if replp is False:
                filter.min_repln(engine, replicate_model, marker.id, 1)
            else:
                filter.min_replp(engine, replicate_model, marker.id, 3)
            filter.indel(False)
            df_codon_stop_per_genetic_code = self.codon_stop_dataframe(genetic_code_tsv)
            filter.codon_stop(df_codon_stop_per_genetic_code, 2, False)