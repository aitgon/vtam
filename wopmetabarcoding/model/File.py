from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer


class File(Base):
    __tablename__ = 'File'

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(150), nullable=False)
    # Todo: (Opt) Create model "Run"?
    run_name = Column(String(20), nullable=False)
    trimmed_status = Column(String, nullable=False)
