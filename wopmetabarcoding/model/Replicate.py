from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer, ForeignKey
from sqlalchemy.orm import validates


class Replicate(Base):
    __tablename__ = 'Replicate'

    replicate_id = Column(Integer, primary_key=True, autoincrement=True) # Todo: Rename to id
    sample_id = Column(Integer, ForeignKey('Biosample.sample_id'))
    marker_id = Column(Integer, ForeignKey('Marker.marker_id'))
    tag_id = Column(Integer, ForeignKey('Tag.tag_id'))
    file_id = Column(Integer, ForeignKey('File.id'))
    replicate_name = Column(String(50), nullable=False) # Todo: rename to "name"

    @validates('replicate_name')
    def validate_replicatename(self, key, replicate_name):
        if '_' in replicate_name:
            replicate_namebis = replicate_name.replace('_', '')
            assert replicate_namebis.isalnum()
        else:
            assert replicate_name.isalnum()
        return replicate_name