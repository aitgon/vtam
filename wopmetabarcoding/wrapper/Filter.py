from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger
from wopmetabarcoding.wrapper.Filter_function import lfn_per_replicate, lfn_per_variant, lfn_per_readcounts,\
    lfn_per_cutoff, min_repln, min_replp, delete_filtered_variants, pcr_error
from sqlalchemy import select
import pandas


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
        # Input file path
        cutoff_file_tsv = self.input_file(Filter.__input_cutoff_file)
        # Input table models
        marker_model = self.input_table(Filter.__input_table_marker)
        variant_model = self.input_table(Filter.__input_table_variant)
        replicate_model = self.input_table(Filter.__input_table_replicate)
        marker_select = select([marker_model.marker_name])
        marker_obj = engine.execute(marker_select)
        for marker in marker_obj:
            failed_variants = []
            file_name = marker.marker_name + '_sample_count.tsv'
            data_frame = pandas.read_csv(file_name, sep='\t')
            result_filter1 = lfn_per_replicate(engine, replicate_model, variant_model, marker.marker_name, data_frame)
            for sample_replicate in result_filter1:
                variant_list = result_filter1.get(sample_replicate)
                failed_variants += variant_list
            result_filter2 = lfn_per_variant(engine, replicate_model, variant_model, marker.marker_name, data_frame, False)
            for variant in result_filter2:
                variant_list = result_filter2.get(variant)
                failed_variants += variant_list
            result_filter3 = lfn_per_readcounts(engine, replicate_model, variant_model, marker.marker_name, 2, data_frame)
            for variant in result_filter3:
                variant_list = result_filter3.get(variant)
                failed_variants += variant_list
            result_filter4 = lfn_per_cutoff(engine, replicate_model, variant_model, marker.marker_name, data_frame, cutoff_file_tsv, False)
            for variant in result_filter4:
                variant_list = result_filter4.get(variant)
                failed_variants += variant_list
            failed_variants = sorted(set(failed_variants))
            print(failed_variants)
            delete_filtered_variants(session, variant_model, failed_variants)
            session.commit()
            pcr_error(engine, replicate_model, variant_model, 'denied', data_frame, marker.marker_name, False)
            # print(failed_variants)
            # result_min_repln = min_repln(engine, variant_model, replicate_model, marker.marker_name, data_frame)
            # print(result_min_repln)
            # result_min_replp = min_replp(engine, variant_model, replicate_model, marker.marker_name, data_frame, 3)
