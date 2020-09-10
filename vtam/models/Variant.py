from Bio.Seq import Seq
from sqlalchemy import Column, String, Integer
from sqlalchemy.orm import validates
from wopmars.Base import Base

# Compatible with both pre- and post Biopython 1.78:
try:
    from Bio.Alphabet import IUPAC  # No longer availa
except ImportError:
    IUPAC = None

class Variant(Base):
    __tablename__ = __qualname__

    id = Column(Integer, primary_key=True, autoincrement=True)
    sequence = Column(String(250), unique=True, nullable=False)

    @validates('sequence')
    def validate_dna(self, key, value):
        """Validate that this string is a DNA sequence and returns the uppder case"""

        if IUPAC:  # Biopython <1.78
            assert str(
                Seq(value, IUPAC.unambiguous_dna).alphabet) == 'IUPACUnambiguousDNA()'
        else:  # Biopython =>1.78
            assert set(value.upper()).issubset({'A', 'C', 'G', 'T'})


        return value.upper()
