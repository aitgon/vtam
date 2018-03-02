from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer, Boolean
from sqlalchemy.orm import validates


class Sample(Base):
	__tablename__ = 'Sample'

	sample_id = Column(Integer, primary_key=True, autoincrement=True)
	sample_name = Column(String(50), nullable=False)
	positive_control = Column(Boolean, default=False)
	negative_control = Column(Boolean, default=False)

	@validates('sample_name')
	def validate_runname(self, key, sample_name):
		if '_' in sample_name:
			sample_namebis = sample_name.replace('_', '')
			assert sample_namebis.isalnum()
		else:
			assert sample_name.isalnum()
		return sample_name