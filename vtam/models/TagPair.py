from wopmars.Base import Base

from sqlalchemy import Column, String, Integer


class TagPair(Base):
    __tablename__ = 'TagPair'

    id = Column(Integer, primary_key=True, autoincrement=True)
    tag_forward = Column(String(100), nullable=True)
    tag_reverse = Column(String(100), nullable=True)



