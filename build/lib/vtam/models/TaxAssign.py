from wopmars.Base import Base
from sqlalchemy import Column, Integer, ForeignKey, String, UniqueConstraint


class TaxAssign(Base):
    __tablename__ = __qualname__
    __table_args__ = (
        UniqueConstraint('variant_id'),
    )

    id = Column(Integer, primary_key=True, autoincrement=True)
    variant_id = Column(
        Integer,
        ForeignKey(
            "Variant.id",
            onupdate="CASCADE",
            ondelete="CASCADE"),
        nullable=False)
    identity = Column(Integer, nullable=True)
    ltg_rank = Column(Integer, nullable=True)
    ltg_tax_id = Column(Integer, nullable=True)
    ltg_tax_name = Column(String(200), nullable=True)
    blast_db = Column(String(50), nullable=True)
