from sqlalchemy import Boolean
from sqlalchemy import UniqueConstraint
from wopmars.framework.database.Base import Base

from sqlalchemy import Column, String, Integer, ForeignKey


class PassedVariant(Base):
    __tablename__ = "PassedVariant"

    variant_id = Column(Integer, ForeignKey("Variant.id"), primary_key=True)
    passed = Column(Boolean, default=True, nullable=False)
    f1_lfn1_per_replicate = Column(Boolean, default=True, nullable=False)
    f2_lfn2_per_variant = Column(Boolean, default=True, nullable=False)
    f3_lfn2_per_replicate_series = Column(Boolean, default=True, nullable=False)
    f4_lfn3_read_count = Column(Boolean, default=True, nullable=False)
    f5_lfn4_per_variant_with_cutoff = Column(Boolean, default=True, nullable=False)
    f6_lfn4_per_replicate_series_with_cutoff = Column(Boolean, default=True, nullable=False)
    f7_min_repln = Column(Boolean, default=True, nullable=False)
    f8_min_replp = Column(Boolean, default=True, nullable=False)
    f9_pcr_error = Column(Boolean, default=True, nullable=False)
    f10_chimera = Column(Boolean, default=True, nullable=False)
    f10_chimera_borderline = Column(Boolean, default=True, nullable=False)
    f11_renkonen = Column(Boolean, default=True, nullable=False)
    f12_indel = Column(Boolean, default=True, nullable=False)
