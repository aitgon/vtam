import os

from wopmetabarcoding.utils.PathFinder import PathFinder
from wopmetabarcoding.utils.constants import tempdir


def insert_table(session, model, obj):
    try:  # checks if exists Phenotype in db
        session.query(model).filter_by(**obj).one()
    except:  # if not, add
        session.add(model(**obj))


def get_or_create(session, model, **kwargs):
    instance = session.query(model).filter_by(**kwargs).first()
    if instance:
        return instance
    else:
        instance = model(**kwargs)
        session.add(instance)
        session.commit()
        return instance

def create_step_tmp_dir(file):
    this_step_tmp_dir = os.path.join(tempdir, os.path.basename(file))
    PathFinder.mkdir_p(this_step_tmp_dir)
    return this_step_tmp_dir