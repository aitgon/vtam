from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer, ForeignKey
from sqlalchemy import UniqueConstraint
from sqlalchemy.orm import validates

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


class SampleInformation(Base):
    __tablename__ = 'SampleInformation'
    __table_args__ = (
        UniqueConstraint('tag_forward', 'tag_reverse', 'marker_name', 'file_name'),
        UniqueConstraint('marker_name', 'sample_name', 'replicate_name')
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    marker_name = Column(String(50), nullable=False)
    primer_forward = Column(String(100), nullable=True)
    primer_reverse = Column(String(100), nullable=True)
    tag_forward = Column(String(100), nullable=True) # Todo: Create Tag model
    tag_reverse = Column(String(100), nullable=True)
    file_name = Column(String(150), nullable=False) # Todo: should be foreign key to File model ID
    run_name = Column(String(20), nullable=False)
    sample_name = Column(String(50), nullable=False)
    replicate_name = Column(String, nullable=False)

    @validates('file_name')
    def validates_filename(self, key, path):
        assert ' ' not in path
        return path

    @validates('marker_name', 'sample_name', 'run_name')
    def validates_names(self, key, name):
        if '_' in name:
            namebis = name.replace('_', '')
            assert namebis.isalnum()
        else:
            assert name.isalnum()
        return name

    @validates('tag_forward', 'primer_forward', 'tag_reverse', 'primer_reverse')
    def validate_sequences(self, key, sequence):
        if sequence != "":
            assert Seq(sequence, IUPAC.ambiguous_dna)
        return sequence



