from sqlalchemy import UniqueConstraint
from wopmars.Base import Base

from sqlalchemy import Column, String, Integer, ForeignKey


class Variant(Base):
    __tablename__ = __qualname__

    id = Column(Integer, primary_key=True, autoincrement=True)
    sequence = Column(String(250), unique=True, nullable=False)
