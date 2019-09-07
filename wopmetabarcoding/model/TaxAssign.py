from wopmars.framework.database.Base import Base

from sqlalchemy import Column, Integer, ForeignKey
from sqlalchemy import UniqueConstraint


class TaxAssign(Base):
    __tablename__ = "TaxAssign"
    __table_args__ = (
        UniqueConstraint('variant_id'),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    variant_id = Column(Integer, ForeignKey("Variant.id", onupdate="CASCADE", ondelete="CASCADE"), nullable=False)
    identity = Column(Integer, nullable=True)
    ltg_rank = Column(Integer, nullable=True)
    ltg_tax_id = Column(Integer, nullable=True)










