from wopmars.Base import Base
from sqlalchemy import Column, String, Integer


class Marker(Base):
    __tablename__ = __qualname__

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(50), nullable=False)
