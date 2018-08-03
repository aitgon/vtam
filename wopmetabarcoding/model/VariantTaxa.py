from sqlalchemy import Boolean
from sqlalchemy import UniqueConstraint
from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer, ForeignKey


class VariantTaxa(Base):
    __tablename__ = "VariantTaxa"

    variant_id = Column(Integer, ForeignKey("Variant.id"), primary_key=True)
    tax_id = Column(Integer, nullable=True)