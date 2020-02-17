from wopmars.Base import Base

from sqlalchemy import Column, String, Integer, Boolean, select


class Biosample(Base):
    __tablename__ = __qualname__

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(50), nullable=False)
    # positive_control = Column(Boolean, default=False, nullable=True)
    # negative_control = Column(Boolean, default=False, nullable=True)
