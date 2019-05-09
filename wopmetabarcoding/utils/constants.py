import os
import tempfile

#####################
#
# Workflow tmp dir
from wopmetabarcoding.utils.PathFinder import PathFinder

VTAMTMPDIR = None
if 'VTAMTMPDIR' in os.environ:
    VTAMTMPDIR = os.environ['VTAMTMPDIR']
    PathFinder.mkdir_p(VTAMTMPDIR)
tempdir = tempfile.mkdtemp(dir=VTAMTMPDIR)
#
#####################

# Old tax_assign parameters
order = [100.0, 97.0, 95.0, 90.0]
# order = [100.0, 97.0, 95.0, 90.0, 85.0, 80.0]
taxonomic_levels = {"family": 5, "order": 4, "genus": 3, "species": 2, "subspecies": 1}

# New tax_assign parameters
rank_hierarchy =['no rank', 'phylum', 'superclass', 'class', 'subclass', 'infraclass', 'superorder', 'order',
                 'suborder', 'infraorder', 'family', 'subfamily', 'genus', 'subgenus', 'species', 'subspecies']


data_dir = os.path.join("{}/tmp/vtam".format(os.environ['HOME']))
wop_dir = os.path.join("{}/Software/repositories/wopmetabarcodin".format(os.environ['HOME']))

public_data_dir = "http://pedagogix-tagc.univ-mrs.fr/~gonzalez/vtam/"

identity_list = [100, 99, 97, 95, 90, 85, 80, 75, 70]
