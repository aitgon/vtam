from wopmars.Base import Base

from sqlalchemy import Column, String, Integer, UniqueConstraint, Boolean, ForeignKey


class Fasta(Base):
    __tablename__ = __qualname__
    __table_args__ = (
        UniqueConstraint('name'),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(150), nullable=False)
    run_id = Column(Integer, ForeignKey("Run.id", onupdate="CASCADE", ondelete="CASCADE"), nullable=False)
    # is_trimmed = Column(Boolean, nullable=False)
