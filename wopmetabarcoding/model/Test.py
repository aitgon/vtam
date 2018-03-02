from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer, ForeignKey
from sqlalchemy.orm import relationship


class Test(Base):
	__tablename__ = 'Test'

	id = Column(Integer, primary_key=True, autoincrement=True)
	primer_id = Column(Integer, ForeignKey('PrimerPair.id'))
	primerpair = relationship("PrimerPair", )
	marker_id = Column(Integer, ForeignKey('Marker.marker_id'))