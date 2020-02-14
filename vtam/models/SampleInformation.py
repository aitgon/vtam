from wopmars.Base import Base

from sqlalchemy import Column, String, Integer, ForeignKey
from sqlalchemy import UniqueConstraint
from sqlalchemy.orm import validates

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


class SampleInformation(Base):
    __tablename__ = __qualname__
    __table_args__ = (
        UniqueConstraint('run_id', 'marker_id', 'biosample_id', 'replicate'),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    run_id = Column(Integer, ForeignKey("Run.id", onupdate="CASCADE", ondelete="CASCADE"), nullable=False)
    marker_id = Column(Integer, ForeignKey("Marker.id", onupdate="CASCADE", ondelete="CASCADE"), nullable=False)
    # primer_forward = Column(String(100), nullable=False)
    # primer_reverse = Column(String(100), nullable=False)
    # tag_forward = Column(String(100), nullable=False)
    # tag_reverse = Column(String(100), nullable=False)
    fasta_id = Column(Integer, ForeignKey("Fasta.id", onupdate="CASCADE", ondelete="CASCADE"), nullable=False)
    biosample_id = Column(Integer, ForeignKey("Biosample.id", onupdate="CASCADE", ondelete="CASCADE"), nullable=False)
    replicate = Column(Integer, nullable=False)

    @validates('name', 'sample_name', 'run_name')
    def validates_names(self, key, name):
        if '_' in name:
            namebis = name.replace('_', '')
            assert namebis.isalnum()
        else:
            assert name.isalnum()
        return name

    # @validates('tag_forward', 'primer_forward', 'tag_reverse', 'primer_reverse')
    # def validate_sequences(self, key, sequence):
    #     if sequence != "":
    #         assert Seq(sequence, IUPAC.ambiguous_dna)
    #     return sequence
