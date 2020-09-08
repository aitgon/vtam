from wopmars.Base import Base
from sqlalchemy import UniqueConstraint, Column, Integer, ForeignKey


class VariantReadCount(Base):
    __tablename__ = __qualname__
    __table_args__ = (
        UniqueConstraint(
            'run_id',
            'variant_id',
            'marker_id',
            'sample_id',
            'replicate'),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    run_id = Column(
        Integer,
        ForeignKey(
            "Run.id",
            onupdate="CASCADE",
            ondelete="CASCADE"),
        nullable=False)
    marker_id = Column(
        Integer,
        ForeignKey(
            "Marker.id",
            onupdate="CASCADE",
            ondelete="CASCADE"),
        nullable=False)
    sample_id = Column(
        Integer,
        ForeignKey(
            "Sample.id",
            onupdate="CASCADE",
            ondelete="CASCADE"),
        nullable=False)
    replicate = Column(Integer, nullable=False)
    variant_id = Column(
        Integer,
        ForeignKey(
            "Variant.id",
            onupdate="CASCADE",
            ondelete="CASCADE"),
        nullable=False)
    read_count = Column(Integer, nullable=False)
