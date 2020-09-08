from wopmars.Base import Base

from sqlalchemy import UniqueConstraint, Boolean, Column, Integer, ForeignKey


class FilterMinReplicateNumber3(Base):
    __tablename__ = __qualname__
    __table_args__ = (
        UniqueConstraint(
            'marker_id',
            'run_id',
            'variant_id',
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
    filter_delete = Column(Boolean, nullable=False)
