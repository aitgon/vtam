from wopmars.Base import Base

from sqlalchemy import Column, String, Integer, ForeignKey
from sqlalchemy.orm import validates


class Replicate(Base):
    __tablename__ = 'Replicate'

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(50), nullable=False)

    @validates('name')
    def validate_replicatename(self, key, replicate_name):
        if '_' in replicate_name:
            replicate_namebis = replicate_name.replace('_', '')
            assert replicate_namebis.isalnum()
        else:
            assert replicate_name.isalnum()
        return replicate_name