from sqlalchemy import UniqueConstraint, Boolean
from wopmars.framework.database.Base import Base

from sqlalchemy import Column, Integer, ForeignKey


class FilterLFN(Base):
    __tablename__ = "FilterLFN"
    __table_args__ = (
        UniqueConstraint('marker_id', 'run_id', 'variant_id', 'biosample_id', 'replicate_id', 'filter_id'),
        # UniqueConstraint('marker_id', 'run_id', 'variant_id', 'biosample_id', 'replicate_id'),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    run_id = Column(Integer, ForeignKey("Run.id", onupdate="CASCADE", ondelete="CASCADE"), nullable=False)
    marker_id = Column(Integer, ForeignKey("Marker.id", onupdate="CASCADE", ondelete="CASCADE"), nullable=False)
    biosample_id = Column(Integer, ForeignKey("Biosample.id", onupdate="CASCADE", ondelete="CASCADE"), nullable=False)
    replicate_id = Column(Integer, ForeignKey("Replicate.id", onupdate="CASCADE", ondelete="CASCADE"), nullable=False)
    variant_id = Column(Integer, ForeignKey("Variant.id", onupdate="CASCADE", ondelete="CASCADE"), nullable=False)
    read_count = Column(Integer, nullable=False)
    filter_id = Column(Integer, nullable=False)
    filter_delete = Column(Boolean, nullable=False)
