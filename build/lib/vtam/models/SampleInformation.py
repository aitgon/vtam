from wopmars.Base import Base

from sqlalchemy import Column, Integer, ForeignKey
from sqlalchemy import UniqueConstraint
from sqlalchemy.orm import validates


class SampleInformation(Base):
    __tablename__ = __qualname__
    __table_args__ = (
        UniqueConstraint(
            'run_id',
            'marker_id',
            'sample_id',
            'replicate',
            'sortedreadfile_id'),
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
    sortedreadfile_id = Column(
        Integer,
        ForeignKey(
            "SortedReadFile.id",
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

    @validates('name', 'sample_name', 'run_name')
    def validates_names(self, key, name):
        if '_' in name:
            namebis = name.replace('_', '')
            assert namebis.isalnum()
        else:
            assert name.isalnum()
        return name
