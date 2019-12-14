from sqlalchemy import UniqueConstraint
from wopmars.Base import Base

from sqlalchemy import Column, Integer, ForeignKey


class VariantReadCount(Base):
    __tablename__ = "VariantReadCount"
    __table_args__ = (
        UniqueConstraint('run_id', 'variant_id', 'marker_id', 'biosample_id', 'replicate'),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    run_id = Column(Integer, ForeignKey("Run.id", onupdate="CASCADE", ondelete="CASCADE"), nullable=False)
    marker_id = Column(Integer, ForeignKey("Marker.id", onupdate="CASCADE", ondelete="CASCADE"), nullable=False)
    biosample_id = Column(Integer, ForeignKey("Biosample.id", onupdate="CASCADE", ondelete="CASCADE"), nullable=False)
    # replicate = Column(Integer, ForeignKey("Replicate.id", onupdate="CASCADE", ondelete="CASCADE"), nullable=False)
    replicate = Column(Integer, nullable=False)
    variant_id = Column(Integer, ForeignKey("Variant.id", onupdate="CASCADE", ondelete="CASCADE"), nullable=False)
    read_count = Column(Integer, nullable=False)
