# -*- coding: utf-8 -*-

import os
import pathlib
import urllib.request

from vtam.utils.constants import public_data_dir


class CommandBlastCOI(object):

    def __init__(self, coi_blast_db_dir, map_taxids_tsv=None):
        self.coi_blast_db = coi_blast_db_dir
        self.map_taxids_tsv = map_taxids_tsv

    ##########################################################
    #
    # Define/create map_taxids and coi_blast_db
    #
    ##########################################################

    def download(self):
        """
        These function is used to define and return the output of the COI Blast database directory.

        If the COI_BLAST_DB environment variable is passed, then that output will be used as COI Blast database directory.
        Otherwise the COI Blast database will be downloaded from the VTAM public data dir to the VTAM local data directory

            Updated:
                May 31, 2019

            Returns:
                String: The output to the taxonomy.sqlite database
        """
        pathlib.Path(os.path.join(self.coi_blast_db)).mkdir(exist_ok=True)

        ####################################################################
        #
        ####################################################################

        basename = 'coi_blast_db.nhr'
        coi_db_url = os.path.join(public_data_dir, basename)
        coi_db_path = os.path.join(self.coi_blast_db, basename)
        if not os.path.isfile(os.path.join(coi_db_path)):
            urllib.request.urlretrieve(coi_db_url, coi_db_path)

        ####################################################################
        #
        ####################################################################

        basename = 'coi_blast_db.nin'
        coi_db_url = os.path.join(public_data_dir, basename)
        coi_db_path = os.path.join(self.coi_blast_db, basename)
        if not os.path.isfile(os.path.join(coi_db_path)):
            urllib.request.urlretrieve(coi_db_url, coi_db_path)

        ####################################################################
        #
        ####################################################################

        basename = 'coi_blast_db.nog'
        coi_db_url = os.path.join(public_data_dir, basename)
        coi_db_path = os.path.join(self.coi_blast_db, basename)
        if not os.path.isfile(os.path.join(coi_db_path)):
            urllib.request.urlretrieve(coi_db_url, coi_db_path)

        ####################################################################
        #
        ####################################################################

        basename = 'coi_blast_db.nsd'
        coi_db_url = os.path.join(public_data_dir, basename)
        coi_db_path = os.path.join(self.coi_blast_db, basename)
        if not os.path.isfile(os.path.join(coi_db_path)):
            urllib.request.urlretrieve(coi_db_url, coi_db_path)

        ####################################################################
        #
        ####################################################################

        basename = 'coi_blast_db.nsi'
        coi_db_url = os.path.join(public_data_dir, basename)
        coi_db_path = os.path.join(self.coi_blast_db, basename)
        if not os.path.isfile(os.path.join(coi_db_path)):
            urllib.request.urlretrieve(coi_db_url, coi_db_path)

        ####################################################################
        #
        ####################################################################

        basename = "coi_blast_db.nsq"
        coi_db_url = os.path.join(public_data_dir, basename)
        coi_db_path = os.path.join(self.coi_blast_db, basename)
        if not os.path.isfile(os.path.join(coi_db_path)):
            urllib.request.urlretrieve(coi_db_url, coi_db_path)

        ####################################################################
        #
        ####################################################################
