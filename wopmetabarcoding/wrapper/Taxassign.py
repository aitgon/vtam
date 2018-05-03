from wopmars.framework.database. tables.ToolWrapper import ToolWrapper
from sqlalchemy import select


class Taxassign(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.Taxassign"
    }
    __input_table_marker = "Marker"
    __input_table_variant = "Variant"
    __input_file_taxassign_db = "taxassign_db_fasta"

    def specify_input_table(self):
        return [
            Taxassign.__input_table_marker,
            Taxassign.__input_table_variant
        ]

    def specify_input_file(self):
        return [
            Taxassign.__input_file_taxassign_db
        ]

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        # Input table models
        marker_model = self.input_table(Taxassign.__input_table_marker)
        variant_model = self.input_table(Taxassign.__input_table_variant)
        # Input files
        taxassign_db_fasta = self.input_file(Taxassign.__input_file_taxassign_db)
        marker_select = select([marker_model.id])
        marker_obj = engine.execute(marker_select)
        for marker in marker_obj:
            variant_select = select([variant_model.variant_id, variant_model.sequence]).where(variant_model.marker_id == marker.id)
            variant_obj = engine.execute(variant_select)
