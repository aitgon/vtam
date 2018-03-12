from wopmars.framework.database.Base import Base

from sqlalchemy.orm import validates
from sqlalchemy import Column, String, Integer


class ReverseTrimmed(Base):
	__tablename__ = "ReverseTrimmed"

	sequence_id = Column(String, primary_key=True)
	target = Column(String, nullable=False)
	tl = Column(Integer, nullable=False)
	qilo = Column(Integer, nullable=False)
	qihi = Column(Integer, nullable=False)
	tilo = Column(Integer, nullable=False)
	tihi = Column(Integer, nullable=False)
	qrow = Column(String, nullable=False)
	read = Column(String, nullable=True)
	reverse = Column(String, nullable=True)
	filename = Column(String, nullable=True)

	@validates('tilo')
	def validate_first_pos(self, key, n):
		assert n == "1"
		return n

	@validates('tihi', 'tl')
	def validate_last_pos(self, key, n):
		tihi_value = ""
		if key == 'tihi':
			tihi_value = n
		if key == "tl":
			tihi_value == n
		return n

	@validates('target', 'qrow')
	def validate_tag_match(self, key, sequence):
		tag_sequence = ""
		if key == 'target':
			for character in sequence:
				if character.islower():
					tag_sequence += character
		if key == 'qrow':
			assert tag_sequence in sequence
		return sequence