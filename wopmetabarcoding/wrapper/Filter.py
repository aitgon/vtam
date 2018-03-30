from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger
from wopmetabarcoding.wrapper.Filter_function import filter1, filter2
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
            file_name = marker.marker_name + '_sample_count.tsv'
            data_frame = pandas.read_csv(file_name, sep='\t')
            result_filter1 = filter2(engine, replicate_model, variant_model, marker.marker_name, file_name, data_frame)
            for element in result_filter1:
                failed = result_filter1.get(element)
                if len(failed) != 0:
                    print(failed)








