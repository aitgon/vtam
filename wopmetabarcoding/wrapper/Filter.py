from wopmars.utils.Logger import Logger
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger

from wopmetabarcoding.utils.constants import tempdir
from wopmetabarcoding.wrapper.FilterUtilities import lfn1_per_replicate, lfn2_per_variant, lfn2_per_replicate_series, \
    lfn3_read_count, lfn4_per_variant_with_cutoff, lfn4_per_replicate_series_with_cutoff, delete_filtered_variants, \
    min_repln, chimera, min_replp, pcr_error, renkoken
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

    def specify_input_table(self):
        return [
            Filter.__input_table_marker,
            Filter.__input_table_variant,
            Filter.__input_table_replicate
        ]

    def specify_input_file(self):
        return[
            Filter.__input_cutoff_file
        ]

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        conn = engine.connect()
        #
        # Input file path
        cutoff_file_tsv = self.input_file(Filter.__input_cutoff_file)
        #
        # Input table models
        marker_model = self.input_table(Filter.__input_table_marker)
        variant_model = self.input_table(Filter.__input_table_variant)
        replicate_model = self.input_table(Filter.__input_table_replicate)
        #
        marker_select = select([marker_model.id, marker_model.name])
        marker_obj = engine.execute(marker_select)
        for marker in marker_obj:
            failed_variants = []
            file_name = os.path.join(tempdir, marker.name + '_sample_count.tsv')
            df = pandas.read_csv(file_name, sep='\t')
            lfn_per_replicate_threshold = 0.001
            lfn_per_variant_threshold = 0.001
            lfn_per_replicate_series_threshold = 0.001
            lfn_read_count_threshold = 2
            #
            Logger.instance().info("Launching LFN filter:")
            failed_variant_list_lfn1_per_replicate = lfn1_per_replicate(df, lfn_per_replicate_threshold)
            lfn_per_replicate_series = False
            if not lfn_per_replicate_series:
                failed_variant_list_lfn2 = lfn2_per_variant(df, lfn_per_variant_threshold)
            else:
                failed_variant_list_lfn2 = lfn2_per_replicate_series(df, lfn_per_replicate_series_threshold)
            failed_variant_list_lfn3_read_count = lfn3_read_count(df, lfn_read_count_threshold)
            if not lfn_per_replicate_series:
                failed_variant_list_lfn4 = lfn4_per_variant_with_cutoff(df, cutoff_file_tsv)
            else:
                failed_variant_list_lfn4 = lfn4_per_replicate_series_with_cutoff(df, cutoff_file_tsv)
            #
            delete_filtered_variants(
                df, failed_variant_list_lfn1_per_replicate, failed_variant_list_lfn2,
                failed_variant_list_lfn3_read_count, failed_variant_list_lfn4
            )
            # min_repln(df, 2)
            Logger.instance().info("Launching Repeatability variants filters:")
            replp = False
            if replp is False:
                min_repln(engine, replicate_model, marker.id, df, 2)
            else:
                min_replp(engine, replicate_model, marker.id, df, 3)
            Logger.instance().info("Launching PCR error filter:")
            pcr_error(engine, replicate_model, variant_model, df, marker.id, 0.5, False)
            # Logger.instance().info("Launching chimera filter:")
            # chimera(engine, replicate_model, variant_model, df, marker.id, 'sample_replicate')
            # if replp is False:
            #     min_repln(engine, variant_model, replicate_model, marker.id, df, 1)
            # else:
            #     min_replp(engine, variant_model, replicate_model, marker.id, df, 3)
            # renkoken(engine, replicate_model, variant_model, df, marker.id)