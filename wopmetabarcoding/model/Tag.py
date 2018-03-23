from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer
from sqlalchemy.orm import validates

import Bio
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import Alphabet

# Todo: Rename to TagPair
# Todo: Maybe will need Tag model with each Tag including fwd and rev
class Tag(Base):
    __tablename__ = 'Tag'

    tag_id = Column(Integer, primary_key=True, autoincrement=True)
    tag_forward = Column(String(100), nullable=True)
    tag_reverse = Column(String(100), nullable=True)



