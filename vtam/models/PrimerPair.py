from wopmars.Base import Base

from sqlalchemy import Column, String, Integer


class PrimerPair(Base):
    __tablename__ = __qualname__

    id = Column(Integer, primary_key=True, autoincrement=True)
    primer_forward = Column(String(100), nullable=True)
    primer_reverse = Column(String(100), nullable=True)