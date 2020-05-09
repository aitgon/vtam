from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from sqlalchemy import Column, String, Integer
from sqlalchemy.orm import validates
from wopmars.Base import Base


class Variant(Base):
    __tablename__ = __qualname__

    id = Column(Integer, primary_key=True, autoincrement=True)
    sequence = Column(String(250), unique=True, nullable=False)

    @validates('sequence')
    def validate_dna(self, key, value):
        """Validate that this string is a DNA sequence and returns the uppder case"""
        assert str(Seq(value, IUPAC.unambiguous_dna).alphabet) == 'IUPACUnambiguousDNA()'
        return value.upper()
