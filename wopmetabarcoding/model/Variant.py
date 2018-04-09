from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer


class Variant(Base):
    __tablename__ = "Variant"

    variant_id = Column(String, primary_key=True)
    # Todo: There is a marker model, so here there should be a marker_id
    # Todo2: According to "data_model_v2.dia" Variant model is composed of: sequence and replicate_id
    marker = Column(String, nullable=False) # Todo: remplace par marker_id
    sequence = Column(String, nullable=False)