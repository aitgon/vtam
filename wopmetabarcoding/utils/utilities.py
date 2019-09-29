import inspect
import os
import tempfile
import urllib

from wopmetabarcoding.utils.PathManager import PathManager
from wopmetabarcoding.utils.constants import public_data_dir
from wopmetabarcoding.utils.Logger import Logger


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
    PathManager.mkdir_p(this_step_tmp_dir)
    return this_step_tmp_dir


##########################################################
#
# Define/create map_taxids and coi_blast_db
#
##########################################################
def download_coi_db():
    """
    These function is used to define and return the path of the COI Blast database directory.

    If the COI_BLAST_DB environment variable is passed, then that path will be used as COI Blast database directory.
    Otherwise the COI Blast database will be downloaded from the VTAM public data dir to the VTAM local data directory

        Updated:
            May 31, 2019

        Args:
            None

        Returns:
            String: The path to the taxonomy.sqlite database
    """
    if 'COI_BLAST_DB' in os.environ:
        return os.environ['COI_BLAST_DB']
    ####
    # vtam_data_dir and coi_blast_db dir
    ####
    vtam_data_dir = create_vtam_data_dir()
    coi_blast_db_dir = os.path.join(vtam_data_dir, 'coi_blast_db')
    PathManager.mkdir_p(os.path.join(coi_blast_db_dir))
    ####
    # map_taxids.tsv
    ####
    map_taxids_tsv_url = os.path.join(public_data_dir, "map_taxids.tsv")
    map_taxids_tsv_path = os.path.join(vtam_data_dir, 'map_taxids.tsv')
    if not os.path.isfile(os.path.join(map_taxids_tsv_path)):
        urllib.request.urlretrieve(map_taxids_tsv_url, map_taxids_tsv_path)
    ####
    # coi_db.nhr
    ####
    coi_db_url = os.path.join(public_data_dir, "coi_blast_db", "coi_db.nhr")
    coi_db_path = os.path.join(vtam_data_dir, 'coi_blast_db', 'coi_db.nhr')
    if not os.path.isfile(os.path.join(coi_db_path)):
        urllib.request.urlretrieve(coi_db_url, coi_db_path)
    ####
    # coi_db.nin
    ####
    coi_db_url = os.path.join(public_data_dir, "coi_blast_db", "coi_db.nin")
    coi_db_path = os.path.join(vtam_data_dir, 'coi_blast_db', 'coi_db.nin')
    if not os.path.isfile(os.path.join(coi_db_path)):
        urllib.request.urlretrieve(coi_db_url, coi_db_path)
    ####
    # coi_db.nog
    ####
    coi_db_url = os.path.join(public_data_dir, "coi_blast_db", "coi_db.nog")
    coi_db_path = os.path.join(vtam_data_dir, 'coi_blast_db', 'coi_db.nog')
    if not os.path.isfile(os.path.join(coi_db_path)):
        urllib.request.urlretrieve(coi_db_url, coi_db_path)
    ####
    # coi_db.nsd
    ####
    coi_db_url = os.path.join(public_data_dir, "coi_blast_db", "coi_db.nsd")
    coi_db_path = os.path.join(vtam_data_dir, 'coi_blast_db', 'coi_db.nsd')
    if not os.path.isfile(os.path.join(coi_db_path)):
        urllib.request.urlretrieve(coi_db_url, coi_db_path)
    ####
    # coi_db.nsi
    ####
    coi_db_url = os.path.join(public_data_dir, "coi_blast_db", "coi_db.nsi")
    coi_db_path = os.path.join(vtam_data_dir, 'coi_blast_db', 'coi_db.nsi')
    if not os.path.isfile(os.path.join(coi_db_path)):
        urllib.request.urlretrieve(coi_db_url, coi_db_path)
    ####
    # coi_db.nsq
    ####
    coi_db_url = os.path.join(public_data_dir, "coi_blast_db", "coi_db.nsq")
    coi_db_path = os.path.join(vtam_data_dir, 'coi_blast_db', 'coi_db.nsq')
    if not os.path.isfile(os.path.join(coi_db_path)):
        urllib.request.urlretrieve(coi_db_url, coi_db_path)
    #
    return map_taxids_tsv_path, coi_blast_db_dir


##########################################################
#
# Define/create taxonomy.sqlite path
#
##########################################################
def download_taxonomy_sqlite():
    """
    These function is used to define and return the path of the taxonomy.sqlite database.

    If the TAXONOMY_SQLITE environment variable is passed, then that path will be used
    Otherwise the taxonomy.sqlite database will be downloaded from the VTAM public data dir to the VTAM local data directory

        Updated:
            May 31, 2019

        Args:
            None

        Returns:
            String: The path to the taxonomy.sqlite database
    """
    logger.debug(
        "file: {}; line: {}; download_taxonomy_sqlite()".format(__file__,
                                                                  inspect.currentframe().f_lineno,))
    if 'TAXONOMY_SQLITE' in os.environ:
        return os.environ['TAXONOMY_SQLITE']
    ####
    # vtam_data_dir and coi_blast_db dir
    ####
    vtam_data_dir = create_vtam_data_dir()
    PathManager.mkdir_p(os.path.join(vtam_data_dir))
    ####
    # taxonomy.sqlite
    ####
    taxonomy_sqlite_url = os.path.join(public_data_dir, "taxonomy.sqlite")
    taxonomy_sqlite_path = os.path.join(vtam_data_dir, 'taxonomy.sqlite')
    if not os.path.isfile(os.path.join(taxonomy_sqlite_path)):
        urllib.request.urlretrieve(taxonomy_sqlite_url, taxonomy_sqlite_path)
    #
    return taxonomy_sqlite_path



##########################################################
#
# Define/create VTAM_TMP_DIR
#
##########################################################
VTAM_TMP_DIR = None
if 'VTAM_TMP_DIR' in os.environ:
    VTAM_TMP_DIR = os.environ['VTAM_TMP_DIR']
    PathManager.mkdir_p(VTAM_TMP_DIR)
tempdir = tempfile.mkdtemp(dir=VTAM_TMP_DIR)


##########################################################
#
# Define/create VTAM_DATA_DIR
#
##########################################################
def create_vtam_data_dir():
    """
    This function creates a VTAM data directory, where all data of VTAM will be downloaded.

    If the VTAM_DATA_DIR environment variable is passed, then that value will be the VTAM data directory.
    Otherwise the PWD environment will be used.
    Currently these data sets will be downloaded

    - The taxonomy sqlite database
    - The COI Blast database
    - The map_taxids.tsv will the link between accession and taxids

        Updated:
            May 31, 2019

        Args:
            None

        Returns:
            String: The path to the VTAM data directory
    """
    if 'VTAM_DATA_DIR' in os.environ:
        VTAM_DATA_DIR = os.environ['VTAM_DATA_DIR']
        PathManager.mkdir_p(VTAM_DATA_DIR)
        return VTAM_DATA_DIR
    else:
        return os.environ['PWD']
