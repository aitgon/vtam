from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger
# from wopmetabarcoding.wrapper.Filter_function import lfn_filter
from sqlalchemy import select


class Filter(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.Filter"
    }

    # Input tables:
    __input_table_marker = "Marker"
    __input_table_variant = "Variant"
    __input_table_replicate = "Replicate"

    def specify_input_table(self):
        return [
            Filter.__input_table_marker,
            Filter.__input_table_variant,
            Filter.__input_table_replicate
        ]

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind

        # Input table models
        marker_model = self.input_table(Filter.__input_table_marker)
        variant_model = self.input_table(Filter.__input_table_variant)
        replicate_model = self.input_table(Filter.__input_table_replicate)

        replicate_select = select([replicate_model.marker_name, replicate_model.name])
        replicate_obj = engine.execute(replicate_select)

        variant_select = select([variant_model.variant_id, variant_model.sequence])
        variant_obj = engine.execute(variant_select)
        i = 0
        j = 0
        file_name = 'MFZR_sample_count.tsv'
        with open(file_name, 'r') as fin:
            for line in fin:
                for replicate in replicate_obj:
                    for variant in variant_obj:
                        if replicate.marker_name in file_name:
                            if replicate.name and variant.sequence in line:
                                i += 1
                            if replicate.name in line:
                                j += 1




