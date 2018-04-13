from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer, ForeignKey


class Variant(Base):
    __tablename__ = "Variant"

    variant_id = Column(String, nullable=False)
    # Todo: There is a marker_id model, so here there should be a marker_id
    # Todo2: According to "data_model_v2.dia" Variant model is composed of: sequence and replicate_id
    marker_id = Column(Integer, ForeignKey("Marker.id"), nullable=False)
    sequence = Column(String, primary_key=True)
    readcount = Column(Integer, nullable=False)