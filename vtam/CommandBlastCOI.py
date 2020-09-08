import argparse
import os
import pathlib
import tarfile
from urllib import request

from vtam.utils.MyProgressBar import MyProgressBar
from vtam.utils.PathManager import PathManager
from vtam.utils.constants import coi_blast_db_gz_url


class CommandBlastCOI(object):

    def __init__(self, blastdbname='coi_blast_db'):

        self.blastdbname = blastdbname

        self.tempdir = PathManager.instance().get_tempdir()
        pathlib.Path(os.path.join(self.tempdir)).mkdir(exist_ok=True, parents=True)

        self.coi_blast_db_gz_url = os.path.dirname(coi_blast_db_gz_url) + '/{}.tar.gz'.format(self.blastdbname)
        self.coi_blast_db_gz_path = os.path.join(self.tempdir, '{}.tar.gz'.format(self.blastdbname))

    def argparse_checker_blast_coi_blastdbname(self):

        try:
            request.urlopen(self.coi_blast_db_gz_url)
            return self.blastdbname
        except :
            raise argparse.ArgumentTypeError(
                "There is not this COI Blast DB name '{}'. Please check available versions (<version>.tar.gz), for instance 'coi_blast_db', here: '{}'".format(self.blastdbname, os.path.dirname(coi_blast_db_gz_url)))

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

        if not os.path.isfile(self.coi_blast_db_gz_path):
            request.urlretrieve(self.coi_blast_db_gz_url, self.coi_blast_db_gz_path, MyProgressBar())

        tar = tarfile.open(self.coi_blast_db_gz_path)
        pathlib.Path(os.path.join(blastdbdir)).mkdir(exist_ok=True, parents=True)
        tar.extractall(blastdbdir)
        tar.close()
