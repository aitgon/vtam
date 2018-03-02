from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer


class Primer(Base):
	__tablename__ = 'Primer'

	primer_id = Column(Integer, primary_key=True, autoincrement=True)
	primer_forward = Column(String(100), nullable=True)
	primer_reverse = Column(String(100), nullable=True)
