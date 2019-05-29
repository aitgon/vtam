import os
import tempfile
import urllib

from wopmetabarcoding.utils.PathFinder import PathFinder

##########################################################
#
# Define VTAM_TMP_DIR
#
##########################################################
VTAM_TMP_DIR = None
if 'VTAM_TMP_DIR' in os.environ:
    VTAM_TMP_DIR = os.environ['VTAM_TMP_DIR']
    PathFinder.mkdir_p(VTAM_TMP_DIR)
tempdir = tempfile.mkdtemp(dir=VTAM_TMP_DIR)

##########################################################
#
# Define/create VTAM_DATA_DIR
#
##########################################################
def create_vtam_data_dir():
    if 'VTAM_DATA_DIR' in os.environ:
        VTAM_DATA_DIR = os.environ['VTAM_DATA_DIR']
        PathFinder.mkdir_p(VTAM_DATA_DIR)
    else:
        VTAM_DATA_DIR = os.environ['PWD']
    return VTAM_DATA_DIR

##########################################################
#
# Define/create COI_BLAST_DB_DIR
#
##########################################################
def create_coi_blast_db():
    ####
    # vtam_data_dir and coi_blast_db dir
    ####
    vtam_data_dir = create_vtam_data_dir()
    coi_blast_db_dir = os.path.join(vtam_data_dir, 'coi_blast_db')
    PathFinder.mkdir_p(os.path.join(coi_blast_db_dir))
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

    return map_taxids_tsv_path, coi_blast_db_dir

#
#####################

# Old tax_assign parameters
order = [100.0, 97.0, 95.0, 90.0]
# order = [100.0, 97.0, 95.0, 90.0, 85.0, 80.0]
taxonomic_levels = {"family": 5, "order": 4, "genus": 3, "species": 2, "subspecies": 1}

# New tax_assign parameters
rank_hierarchy =['no rank', 'phylum', 'superclass', 'class', 'subclass', 'infraclass', 'superorder', 'order',
                 'suborder', 'infraorder', 'family', 'subfamily', 'genus', 'subgenus', 'species', 'subspecies']
rank_hierarchy_otu_table =['phylum', 'class', 'order', 'family', 'genus', 'species']


wop_dir = os.path.join("{}/Software/repositories/wopmetabarcodin".format(os.environ['HOME']))

public_data_dir = "http://pedagogix-tagc.univ-mrs.fr/~gonzalez/vtam/"

identity_list = [100, 99, 97, 95, 90, 85, 80, 75, 70]
