from wopmars.Base import Base
from sqlalchemy import UniqueConstraint, Float, Column, Integer, ForeignKey


class ReadCountAverageOverReplicates(Base):
    __tablename__ = __qualname__
    __table_args__ = (
        UniqueConstraint('marker_id', 'run_id', 'variant_id', 'sample_id'),
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
    variant_id = Column(
        Integer,
        ForeignKey(
            "Variant.id",
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
    read_count = Column(Integer, nullable=False)
    replicate_count = Column(Integer, nullable=False)
    read_count_average = Column(Float, nullable=False)
