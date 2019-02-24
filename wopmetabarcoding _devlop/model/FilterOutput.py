from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer


class FilterOutput(Base):
    __tablename__ = 'FilterOutput'  # current file

    marker_name = Column(String, primary_key=True)
    dataframe_path = Column(String, nullable=False)
    output_fasta_path = Column(String, nullable=False)