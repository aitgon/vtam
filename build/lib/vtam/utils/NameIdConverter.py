import sys

import sqlalchemy

from vtam.models.FilterChimeraBorderline import FilterChimeraBorderline
from vtam.models.Variant import Variant
from vtam.utils.Logger import Logger


class NameIdConverter:
    """Takes a list of names or IDs and returns the complementeary"""

    def __init__(self, id_name_or_sequence_list, engine):

        self.id_name_or_sequence_list = id_name_or_sequence_list
        self.engine = engine

    def to_ids(self, declarative_model):

        id_lst = []
        with self.engine.connect() as conn:
            for name in self.id_name_or_sequence_list:
                result = conn.execute(
                    sqlalchemy.select([declarative_model.__table__.c.id]).where(
                        declarative_model.__table__.c.name == name)).first()
                if result is None:
                    Logger.instance().error("Name {} not found in table {}".format(name, str(declarative_model.__table__)))
                    sys.exit(1)
                id_lst.append(result[0])
        return id_lst

    def to_names(self, declarative_model):

        nameid_lst = []
        with self.engine.connect() as conn:
            for idx in self.id_name_or_sequence_list:
                result = conn.execute(
                    sqlalchemy.select([declarative_model.__table__.c.name]).where(
                        declarative_model.__table__.c.id == idx)).first()
                if result is None:
                    Logger.instance().error("Id {} not found in table {}".format(idx, str(declarative_model.__table__)))
                    sys.exit(1)
                nameid_lst.append(result[0])
        return nameid_lst

    def variant_id_to_sequence(self):

        sequence_lst = []
        with self.engine.connect() as conn:
            for variant_id in self.id_name_or_sequence_list:
                result = conn.execute(
                    sqlalchemy.select([Variant.__table__.c.sequence]).where(
                        Variant.__table__.c.id == variant_id)).first()
                if result is None:
                    Logger.instance().error("Variant ID {} not found in table {}".format(variant_id, str(Variant.__table__)))
                    sys.exit(1)
                sequence_lst.append(result[0])
        return sequence_lst

    def variant_sequence_to_id(self):

        variant_id_lst = []
        with self.engine.connect() as conn:
            for sequence in self.id_name_or_sequence_list:
                result = conn.execute(
                    sqlalchemy.select([Variant.__table__.c.id]).where(
                        Variant.__table__.c.sequence == sequence)).first()
                if result is None:
                    Logger.instance().error("Sequence {} not found in table {}".format(sequence, str(Variant.__table__)))
                    sys.exit(1)
                variant_id_lst.append(result[0])
        return variant_id_lst

    def variant_id_is_chimera_borderline(self):

        chimera_borderline_lst = []
        with self.engine.connect() as conn:
            for variant_id in self.id_name_or_sequence_list:
                result = conn.execute(
                    sqlalchemy.select([FilterChimeraBorderline.__table__.c.filter_delete]).where(
                        FilterChimeraBorderline.__table__.c.variant_id == variant_id).distinct()).first()
                if result is None:
                    Logger.instance().error("Variant ID {} not found in table FilterChimeraBorderline".format(variant_id))
                    sys.exit(1)
                chimera_borderline_lst.append(result[0])
        return chimera_borderline_lst
