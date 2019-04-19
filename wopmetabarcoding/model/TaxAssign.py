from sqlalchemy import UniqueConstraint
from wopmars.framework.database.Base import Base

from sqlalchemy import Column, Integer, ForeignKey


class TaxAssign(Base):
    __tablename__ = "TaxAssign"
    __table_args__ = (
        UniqueConstraint('variant_id'),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    variant_id = Column(Integer, ForeignKey("Variant.id"), nullable=False)
    tax_id = Column(Integer, nullable=True)
