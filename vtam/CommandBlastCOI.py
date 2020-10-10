import argparse
import os
import pathlib
import tarfile
import urllib.request

from urllib import request
from vtam.utils.PathManager import PathManager
from vtam.utils.constants import coi_blast_db_gz_url1, coi_blast_db_gz_url2, coi_blast_db_gz_url3
from vtam.utils.MyProgressBar import MyProgressBar

class CommandBlastCOI(object):

    def __init__(self, blastdbname='coi_blast_db'):

        self.blastdbname = blastdbname

        self.tempdir = PathManager.instance().get_tempdir()
        pathlib.Path(os.path.join(self.tempdir)).mkdir(exist_ok=True, parents=True)

        coi_blast_db_version_gz_url_dir1 = os.path.dirname(coi_blast_db_gz_url1)
        self.coi_blast_db_gz_url1 = os.path.join(coi_blast_db_version_gz_url_dir1, "{}.tar.gz".format(blastdbname))
        coi_blast_db_gz_url_dir2 = os.path.dirname(coi_blast_db_gz_url2)
        self.coi_blast_db_gz_url2 = os.path.join(coi_blast_db_gz_url_dir2, "{}.tar.gz".format(blastdbname))
        self.coi_blast_db_gz_url_dir3 = os.path.dirname(coi_blast_db_gz_url3)
        self.coi_blast_db_gz_url3 = os.path.join(self.coi_blast_db_gz_url_dir3, "{}.tar.gz".format(blastdbname))

        self.coi_blast_db_gz_path = os.path.join(self.tempdir, '{}.tar.gz'.format(self.blastdbname))

    def argparse_checker_blast_coi_blastdbname(self):

        """Verifies whether this blastdbname exists"""

        try:
            request.urlopen(self.coi_blast_db_gz_url1)
            return self.blastdbname
        except:
            try:
                request.urlopen(self.coi_blast_db_gz_url2)
                return self.blastdbname
            except:
                try:
                    request.urlopen(self.coi_blast_db_gz_url3)
                    return self.blastdbname
                except:
                    raise argparse.ArgumentTypeError(
                        "This is not a valid --blastdbname {}. Valid  values are here : '{}'"
                            .format(self.blastdbname, os.path.dirname(self.coi_blast_db_gz_url_dir3)))

    def download(self, blastdbdir):
        """
        These function is used to define and return the output of the COI Blast database directory.

        If the COI_BLAST_DB environment variable is passed, then that output will be used as COI Blast database directory.
        Otherwise the COI Blast database will be downloaded from the VTAM public data dir to the VTAM local data directory

            Updated:
                Jun 03, 2020

            Returns:
                String: The output to the taxonomy.sqlite database
        """

        # Test first in local dir, otherwise in the remote URLs
        if not os.path.isfile(self.coi_blast_db_gz_path) or pathlib.Path(self.coi_blast_db_gz_path).stat().st_size < 1000000:
            try:
                urllib.request.urlretrieve(self.coi_blast_db_gz_url1, self.coi_blast_db_gz_path, MyProgressBar())
            except Exception:
                try:
                    urllib.request.urlretrieve(self.coi_blast_db_gz_url2, self.coi_blast_db_gz_path, MyProgressBar())
                except Exception:
                    urllib.request.urlretrieve(self.coi_blast_db_gz_url3, self.coi_blast_db_gz_path, MyProgressBar())
        tar = tarfile.open(self.coi_blast_db_gz_path)
        pathlib.Path(os.path.join(blastdbdir)).mkdir(exist_ok=True, parents=True)
        tar.extractall(blastdbdir)
        tar.close()
