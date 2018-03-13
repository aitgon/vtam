from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer


class ObiFasta(Base):
	__tablename__ = "ObiFasta"

	variant_id = Column(String, primary_key=True)
	sample_count = Column(String, nullable=False)
	read_count = Column(Integer, nullable=False)
	read = Column(String, nullable=False)