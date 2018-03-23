from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer, Boolean
from sqlalchemy.orm import validates

# Todo: (Opt) Rename to biosample
class Biosample(Base):
    __tablename__ = 'Biosample'

    id = Column(Integer, primary_key=True, autoincrement=True) # Todo: Rename to 'id'
    name = Column(String(50), nullable=False) # Todo: Rename to 'name'
    positive_control = Column(Boolean, default=False)
    negative_control = Column(Boolean, default=False)