from sqlalchemy import UniqueConstraint, Boolean, Float
from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer, ForeignKey


class FilterConsensus(Base):
    __tablename__ = "FilterConsensus"
    __table_args__ = (
        UniqueConstraint('marker_id', 'run_id', 'variant_id', 'biosample_id'),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    run_id = Column(Integer, ForeignKey("Run.id"), nullable=False)
    marker_id = Column(Integer, ForeignKey("Marker.id"), nullable=False)
    variant_id = Column(Integer, ForeignKey("Variant.id"), nullable=False)
    biosample_id = Column(Integer, ForeignKey("Biosample.id"), nullable=False)
    read_count = Column(Integer, nullable=False)
    replicate_count = Column(Integer, nullable=False)
    read_count_average = Column(Float, nullable=False)
