from sqlalchemy import UniqueConstraint
from wopmars.framework.database.Base import Base

from sqlalchemy import Column, Integer, ForeignKey


class TaxAssign(Base):
    __tablename__ = "Ltg"
    __table_args__ = (
        UniqueConstraint('variant_id','identity'),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    variant_id = Column(Integer, ForeignKey("Variant.id"), nullable=False)
    identity = Column(Integer, ForeignKey("identity"), nullable=False)
    ltg_tax_id = Column(float)
    ltg_rank = Column(str)
    ltg_lineage = Column(float)


