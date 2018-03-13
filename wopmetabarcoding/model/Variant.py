from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer


class Variant(Base):
	__tablename__ = "Variant"

	variant_id = Column(String, primary_key=True)
	marker = Column(String, nullable=False)
	sequence = Column(String, nullable=False)