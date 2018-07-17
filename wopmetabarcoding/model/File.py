from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer, UniqueConstraint, Boolean


class File(Base):
    __tablename__ = 'File'
    __table_args__ = (
        UniqueConstraint('name'),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(150), nullable=False)
    # Todo: (Opt) Create model "Run"?
    run_name = Column(String(20), nullable=False)
    is_trimmed = Column(Boolean, nullable=False)