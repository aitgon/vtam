from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer
from sqlalchemy.orm import validates


class File(Base):
	__tablename__ = 'File'

	file_id = Column(Integer, primary_key=True, autoincrement=True)
	file_name = Column(String(150), nullable=False)
	run_name = Column(String(20), nullable=False)