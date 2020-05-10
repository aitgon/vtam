# -*- coding: utf-8 -*-
import gzip
import os
import pathlib
import tarfile
import urllib.request

from vtam.utils.PathManager import PathManager
from vtam.utils.constants import coi_blast_db_gz_url


class CommandBlastCOI(object):

    def __init__(self, coi_blast_db_dir):
        self.coi_blast_db_dir = coi_blast_db_dir
        self.tempdir = PathManager.instance().get_tempdir()
        pathlib.Path(os.path.join(self.tempdir)).mkdir(exist_ok=True, parents=True)

    ##########################################################
    #
    # Define/create map_taxids and coi_blast_db_dir
    #
    ##########################################################

    def download(self):
        """
        These function is used to define and return the output of the COI Blast database directory.

        If the COI_BLAST_DB environment variable is passed, then that output will be used as COI Blast database directory.
        Otherwise the COI Blast database will be downloaded from the VTAM public data dir to the VTAM local data directory

            Updated:
                Apr 25, 2020

            Returns:
                String: The output to the taxonomy.sqlite database
        """

        coi_blast_db_gz_path = os.path.join(self.tempdir, "coi_blast_db.tar.gz")

        nhr_path = os.path.join(self.coi_blast_db_dir, "coi_blast_db_dir.nhr")
        nin_path = os.path.join(self.coi_blast_db_dir, "coi_blast_db_dir.nin")
        nog_path = os.path.join(self.coi_blast_db_dir, "coi_blast_db_dir.nog")
        nsd_path = os.path.join(self.coi_blast_db_dir, "coi_blast_db_dir.nsd")
        nsi_path = os.path.join(self.coi_blast_db_dir, "coi_blast_db_dir.nsi")
        nsq_path = os.path.join(self.coi_blast_db_dir, "coi_blast_db_dir.nsq")

        if not os.path.isfile(nhr_path) \
                or not os.path.isfile(nin_path) \
                or not os.path.isfile(nog_path) \
                or not os.path.isfile(nsd_path) \
                or not os.path.isfile(nsi_path) \
                or not os.path.isfile(nsq_path):
            urllib.request.urlretrieve(coi_blast_db_gz_url, coi_blast_db_gz_path)

            pathlib.Path(os.path.join(self.coi_blast_db_dir)).mkdir(exist_ok=True, parents=True)

            tar = tarfile.open(coi_blast_db_gz_path)
            tar.extractall(self.coi_blast_db_dir)
            tar.close()

        #     with gzip.open(coi_blast_db_gz_path, 'rb') as fin:
        #         with open(self.taxonomy_tsv_path, 'wb') as fout:
        #             shutil.copyfileobj(fin, fout)
        #     try:
        #         pathlib.Path(taxonomy_tsv_gz_path).unlink()
        #     except FileNotFoundError:
        #         pass
        #
        # ####################################################################
        # #
        # ####################################################################
        #
        # basename = 'coi_blast_db_dir.nhr'
        # coi_db_url = os.tsv_path.join(public_data_dir, basename)
        # coi_db_path = os.tsv_path.join(self.coi_blast_db_dir, basename)
        # if not os.tsv_path.isfile(os.tsv_path.join(coi_db_path)):
        #     urllib.request.urlretrieve(coi_db_url, coi_db_path)
        #
        # ####################################################################
        # #
        # ####################################################################
        #
        # basename = 'coi_blast_db_dir.nin'
        # coi_db_url = os.tsv_path.join(public_data_dir, basename)
        # coi_db_path = os.tsv_path.join(self.coi_blast_db_dir, basename)
        # if not os.tsv_path.isfile(os.tsv_path.join(coi_db_path)):
        #     urllib.request.urlretrieve(coi_db_url, coi_db_path)
        #
        # ####################################################################
        # #
        # ####################################################################
        #
        # basename = 'coi_blast_db_dir.nog'
        # coi_db_url = os.tsv_path.join(public_data_dir, basename)
        # coi_db_path = os.tsv_path.join(self.coi_blast_db_dir, basename)
        # if not os.tsv_path.isfile(os.tsv_path.join(coi_db_path)):
        #     urllib.request.urlretrieve(coi_db_url, coi_db_path)
        #
        # ####################################################################
        # #
        # ####################################################################
        #
        # basename = 'coi_blast_db_dir.nsd'
        # coi_db_url = os.tsv_path.join(public_data_dir, basename)
        # coi_db_path = os.tsv_path.join(self.coi_blast_db_dir, basename)
        # if not os.tsv_path.isfile(os.tsv_path.join(coi_db_path)):
        #     urllib.request.urlretrieve(coi_db_url, coi_db_path)
        #
        # ####################################################################
        # #
        # ####################################################################
        #
        # basename = 'coi_blast_db_dir.nsi'
        # coi_db_url = os.tsv_path.join(public_data_dir, basename)
        # coi_db_path = os.tsv_path.join(self.coi_blast_db_dir, basename)
        # if not os.tsv_path.isfile(os.tsv_path.join(coi_db_path)):
        #     urllib.request.urlretrieve(coi_db_url, coi_db_path)
        #
        # ####################################################################
        # #
        # ####################################################################
        #
        # basename = "coi_blast_db_dir.nsq"
        # coi_db_url = os.tsv_path.join(public_data_dir, basename)
        # coi_db_path = os.tsv_path.join(self.coi_blast_db_dir, basename)
        # if not os.tsv_path.isfile(os.tsv_path.join(coi_db_path)):
        #     urllib.request.urlretrieve(coi_db_url, coi_db_path)
        #
        # ####################################################################
        # #
        # ####################################################################
