from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger
from wopmetabarcoding.wrapper.Filter_function import lfn_per_replicate_filter1#, lfn_per_variant_filter2,\
    # empirical_filter3, lfn_per_variant_sc_filter4
from sqlalchemy import select


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
            file_name = marker.marker_name + '_sample_count.tsv'
            failed_variants_1 = lfn_per_replicate_filter1(
                engine, replicate_model, variant_model, marker.marker_name, file_name
            )
            for element in failed_variants_1:
                failed_list = failed_variants_1.get(element)
                print(element)
                print(failed_list)
            # print('Second step')
            # failed_variants_2 = lfn_per_variant_filter2(
            #     engine, replicate_model, variant_model, marker.marker_name, file_name
            # )
            # for element in failed_variants_2:
            #     failed_list = failed_variants_2.get(element)
            #     print(failed_list)
            # failed_variants_3 = empirical_filter3(
            #     engine, replicate_model, variant_model, marker.marker_name, file_name, 1
            # )
            # failed_variants_4 =lfn_per_variant_sc_filter4(
            #     engine, replicate_model, variant_model, marker.marker_name, file_name, cutoff_file_tsv
            # )







