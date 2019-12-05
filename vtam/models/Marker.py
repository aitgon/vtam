from wopmars.Base import Base

from sqlalchemy import Column, String, Integer
from sqlalchemy.orm import validates

import Bio
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import Alphabet


class Marker(Base):
    __tablename__ = 'Marker'  # current file

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(50), nullable=False)




