from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger
from wopmetabarcoding.wrapper.FilterUtilities import lfn_per_replicate, lfn_per_variant, lfn_per_readcounts,\
    lfn_per_cutoff, min_repln, min_replp, delete_filtered_variants, pcr_error, chimera
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
            file_name = os.path.join('data/output/Sort_reads/', marker.name + '_sample_count.tsv')
            data_frame = pandas.read_csv(file_name, sep='\t')
            result_filter1 = lfn_per_replicate(engine, replicate_model, variant_model, marker.id, data_frame)
            result_filter2 = lfn_per_variant(engine, replicate_model, variant_model, marker.id, data_frame, False)
            result_filter3 = lfn_per_readcounts(engine, replicate_model, variant_model, marker.id, 2, data_frame)
            result_filter4 = lfn_per_cutoff(engine, replicate_model, variant_model, marker.id, data_frame, cutoff_file_tsv, False)
            delete_filtered_variants(
                engine, replicate_model, marker.id, data_frame, result_filter1, result_filter2, result_filter3, result_filter4
            )
            replp = True
            if replp is False:
                min_repln(engine, variant_model, replicate_model, marker.id, data_frame, 1)
            else:
                min_replp(engine, variant_model, replicate_model, marker.id, data_frame, 3)
            chimera(engine, replicate_model, variant_model, data_frame, marker.id, 'sample_replicate')