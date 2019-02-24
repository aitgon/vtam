from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer, UniqueConstraint, Boolean

# TODO rename Fasta to Fasta
class Fasta(Base):
    __tablename__ = 'Fasta'
    __table_args__ = (
        UniqueConstraint('name'),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(150), nullable=False)
    # Todo: (Opt) Create model "Run"?
    run_name = Column(String(20), nullable=False)
    is_trimmed = Column(Boolean, nullable=False)