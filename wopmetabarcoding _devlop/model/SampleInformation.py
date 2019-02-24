from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer, ForeignKey
from sqlalchemy import UniqueConstraint
from sqlalchemy.orm import validates

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


class SampleInformation(Base):
    __tablename__ = 'SampleInformation'
    __table_args__ = (
        UniqueConstraint('tag_forward', 'tag_reverse', 'marker_id', 'fasta_id'),
        UniqueConstraint('marker_id', 'biosample_id', 'replicate_id')
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    marker_id = Column(Integer, ForeignKey("Marker.id"), nullable=False)
    primer_forward = Column(String(100), nullable=False)
    primer_reverse = Column(String(100), nullable=False)
    tag_forward = Column(String(100), nullable=False) # Todo: Create Tag model
    tag_reverse = Column(String(100), nullable=False)
    fasta_id = Column(Integer, ForeignKey("Fasta.id"), nullable=False) # Todo: should be foreign key to Fasta model ID
    run_name = Column(String(20), nullable=False)
    biosample_id = Column(Integer, ForeignKey("Biosample.id"), nullable=False)
    replicate_id = Column(Integer, ForeignKey("Replicate.id"), nullable=False)

    # @validates('file_name')
    # def validates_filename(self, key, path):
    #     assert ' ' not in path
    #     return path

    @validates('name', 'sample_name', 'run_name')
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



