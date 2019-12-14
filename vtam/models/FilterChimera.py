from sqlalchemy import UniqueConstraint, Boolean
from wopmars.Base import Base

from sqlalchemy import Column, String, Integer, ForeignKey


class FilterChimera(Base):
    __tablename__ = "FilterChimera"
    __table_args__ = (
        UniqueConstraint('marker_id', 'run_id', 'variant_id', 'biosample_id', 'replicate'),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    run_id = Column(Integer, ForeignKey("Run.id", onupdate="CASCADE", ondelete="CASCADE"), nullable=False)
    marker_id = Column(Integer, ForeignKey("Marker.id", onupdate="CASCADE", ondelete="CASCADE"), nullable=False)
    biosample_id = Column(Integer, ForeignKey("Biosample.id", onupdate="CASCADE", ondelete="CASCADE"), nullable=False)
    # replicate = Column(Integer, ForeignKey("Replicate.id", onupdate="CASCADE", ondelete="CASCADE"), nullable=False)
    replicate = Column(Integer, nullable=False)
    variant_id = Column(Integer, ForeignKey("Variant.id"), nullable=False)
    read_count = Column(Integer, nullable=False)
    filter_delete = Column(Boolean, nullable=False)
