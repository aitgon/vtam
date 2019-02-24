from sqlalchemy import UniqueConstraint
from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer, ForeignKey


class Variant(Base):
    __tablename__ = "Variant"

    id = Column(Integer, primary_key=True, autoincrement=True)
    sequence = Column(String(250), unique=True, nullable=False)