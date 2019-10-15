import os
import urllib

from vtam.utils.PathManager import PathManager
from vtam.utils.constants import public_data_dir
from vtam.utils.Logger import Logger


def get_or_create(session, model, **kwargs):
    instance = session.query(model).filter_by(**kwargs).first()
    if instance:
        return instance
    else:
        instance = model(**kwargs)
        session.add(instance)
        session.commit()
        return instance


##########################################################
#
# Define/create map_taxids and coi_blast_db
#
##########################################################
def download_coi_db():
    """
    These function is used to define and return the output of the COI Blast database directory.

    If the COI_BLAST_DB environment variable is passed, then that output will be used as COI Blast database directory.
    Otherwise the COI Blast database will be downloaded from the VTAM public data dir to the VTAM local data directory

        Updated:
            May 31, 2019

        Args:
            None

        Returns:
            String: The output to the taxonomy.sqlite database
    """
    if 'COI_BLAST_DB' in os.environ:
        return os.environ['COI_BLAST_DB']
    ####
    # vtam_data_dir and coi_blast_db dir
    ####
    tempdir = PathManager.instance().get_tempdir()
    coi_blast_db_dir = os.path.join(tempdir, 'coi_blast_db')
    PathManager.mkdir_p(os.path.join(coi_blast_db_dir))
    ####
    # map_taxids.tsv
    ####
    map_taxids_tsv_url = os.path.join(public_data_dir, "map_taxids.tsv")
    map_taxids_tsv_path = os.path.join(tempdir, 'map_taxids.tsv')
    if not os.path.isfile(os.path.join(map_taxids_tsv_path)):
        urllib.request.urlretrieve(map_taxids_tsv_url, map_taxids_tsv_path)
    ####
    # coi_db.nhr
    ####
    coi_db_url = os.path.join(public_data_dir, "coi_blast_db", "coi_db.nhr")
    coi_db_path = os.path.join(tempdir, 'coi_blast_db', 'coi_db.nhr')
    if not os.path.isfile(os.path.join(coi_db_path)):
        urllib.request.urlretrieve(coi_db_url, coi_db_path)
    ####
    # coi_db.nin
    ####
    coi_db_url = os.path.join(public_data_dir, "coi_blast_db", "coi_db.nin")
    coi_db_path = os.path.join(tempdir, 'coi_blast_db', 'coi_db.nin')
    if not os.path.isfile(os.path.join(coi_db_path)):
        urllib.request.urlretrieve(coi_db_url, coi_db_path)
    ####
    # coi_db.nog
    ####
    coi_db_url = os.path.join(public_data_dir, "coi_blast_db", "coi_db.nog")
    coi_db_path = os.path.join(tempdir, 'coi_blast_db', 'coi_db.nog')
    if not os.path.isfile(os.path.join(coi_db_path)):
        urllib.request.urlretrieve(coi_db_url, coi_db_path)
    ####
    # coi_db.nsd
    ####
    coi_db_url = os.path.join(public_data_dir, "coi_blast_db", "coi_db.nsd")
    coi_db_path = os.path.join(tempdir, 'coi_blast_db', 'coi_db.nsd')
    if not os.path.isfile(os.path.join(coi_db_path)):
        urllib.request.urlretrieve(coi_db_url, coi_db_path)
    ####
    # coi_db.nsi
    ####
    coi_db_url = os.path.join(public_data_dir, "coi_blast_db", "coi_db.nsi")
    coi_db_path = os.path.join(tempdir, 'coi_blast_db', 'coi_db.nsi')
    if not os.path.isfile(os.path.join(coi_db_path)):
        urllib.request.urlretrieve(coi_db_url, coi_db_path)
    ####
    # coi_db.nsq
    ####
    coi_db_url = os.path.join(public_data_dir, "coi_blast_db", "coi_db.nsq")
    coi_db_path = os.path.join(tempdir, 'coi_blast_db', 'coi_db.nsq')
    if not os.path.isfile(os.path.join(coi_db_path)):
        urllib.request.urlretrieve(coi_db_url, coi_db_path)
    #
    return map_taxids_tsv_path, coi_blast_db_dir


##########################################################
#
# Convert DF to list of dictionaries to use in an sqlalchemy core insert
#
##########################################################
def filter_delete_df_to_dict(filter_df):
    """Convert DF to list of dictionaries to use in an sqlalchemy core insert"""
    records = []
    for row in filter_df.itertuples():
        run_id = row.run_id
        marker_id = row.marker_id
        biosample_id = row.biosample_id
        replicate_id = row.replicate_id
        variant_id = row.variant_id
        read_count = row.read_count
        filter_delete = row.filter_delete
        if not 'filter_id' in dir(row): # Default not filter
            records.append({'run_id': run_id, 'marker_id': marker_id,
                                                     'variant_id': variant_id, 'biosample_id': biosample_id,
                                                     'replicate_id': replicate_id, 'read_count': read_count,
                                                     'filter_delete': filter_delete})
        else: # Only used for filterlfn where fil
            filter_id = row.filter_id
            records.append({'run_id': run_id, 'marker_id': marker_id,
                                                     'variant_id': variant_id, 'biosample_id': biosample_id,
                                                     'replicate_id': replicate_id, 'read_count': read_count,
                                                     'filter_id': filter_id, 'filter_delete': filter_delete})
    return records
