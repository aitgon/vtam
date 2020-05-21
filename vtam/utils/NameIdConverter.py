import sqlalchemy

from vtam.models.FilterChimeraBorderline import FilterChimeraBorderline
from vtam.models.Variant import Variant


class NameIdConverter:
    """Takes a list of names or IDs and returns the complementeary"""

    def __init__(self, id_name_or_sequence_list, engine):

        self.id_name_or_sequence_list = id_name_or_sequence_list
        self.engine = engine

    def to_ids(self, declarative_model):

        with self.engine.connect() as conn:
            id_list = list(map(lambda x: conn.execute(
                sqlalchemy.select([declarative_model.__table__.c.id]).where(
                    declarative_model.__table__.c.name == x)).first()[0], self.id_name_or_sequence_list))
        return id_list

    def to_names(self, declarative_model):

        with self.engine.connect() as conn:
            id_list = list(map(lambda x: conn.execute(
                sqlalchemy.select([declarative_model.__table__.c.name]).where(
                    declarative_model.__table__.c.id == x)).first()[0], self.id_name_or_sequence_list))
        return id_list

    def variant_id_to_sequence(self):

        with self.engine.connect() as conn:
            id_list = list(map(lambda x: conn.execute(
                sqlalchemy.select([Variant.__table__.c.sequence]).where(
                    Variant.__table__.c.id == x)).first()[0], self.id_name_or_sequence_list))
        return id_list

    def variant_sequence_to_id(self):

        with self.engine.connect() as conn:
            id_list = list(map(lambda x: conn.execute(
                sqlalchemy.select([Variant.__table__.c.id]).where(
                    Variant.__table__.c.sequence == x)).first()[0], self.id_name_or_sequence_list))
        return id_list

    def variant_id_is_chimera_borderline(self):

        with self.engine.connect() as conn:
            # import pdb; pdb.set_trace()
            id_list = list(map(lambda x: conn.execute(
                sqlalchemy.select([FilterChimeraBorderline.__table__.c.filter_delete]).where(
                    FilterChimeraBorderline.__table__.c.variant_id == x).distinct()).first()[0], self.id_name_or_sequence_list))
        return id_list
