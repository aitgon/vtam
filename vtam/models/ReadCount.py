from wopmars.Base import Base
from sqlalchemy import Column, String, Integer


class ReadCount(Base):
    __tablename__ = __qualname__

    sequence = Column(String, primary_key=True)
    count = Column(Integer, nullable=False)
