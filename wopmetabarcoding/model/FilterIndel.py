from sqlalchemy import UniqueConstraint, Boolean
from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer, ForeignKey


class FilterIndel(Base):
    __tablename__ = "FilterIndel"
    __table_args__ = (
        UniqueConstraint('marker_id', 'run_id', 'variant_id', 'biosample_id', 'replicate_id'),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    run_id = Column(Integer, ForeignKey("Run.id"), nullable=False)
    marker_id = Column(Integer, ForeignKey("Marker.id"), nullable=False)
    variant_id = Column(Integer, ForeignKey("Variant.id"), nullable=False)
    biosample_id = Column(Integer, ForeignKey("Biosample.id"), nullable=False)
    replicate_id = Column(Integer, ForeignKey("Replicate.id"), nullable=False)
    read_count = Column(Integer, nullable=False)
    filter_delete = Column(Boolean, nullable=False)