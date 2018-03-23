from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer
from sqlalchemy.orm import validates

import Bio
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import Alphabet

# Todo: Rename to TagPair
# Todo: Maybe will need TagPair model with each TagPair including fwd and rev
class TagPair(Base):
    __tablename__ = 'TagPair'

    id = Column(Integer, primary_key=True, autoincrement=True)
    tag_forward = Column(String(100), nullable=True)
    tag_reverse = Column(String(100), nullable=True)



