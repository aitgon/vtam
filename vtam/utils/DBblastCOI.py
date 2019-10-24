#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import urllib.request

from vtam.utils.PathManager import PathManager
from vtam.utils.constants import public_data_dir

class DBblastCOI(object):

    def __init__(self, coi_blast_db, map_taxids_tsv=None):
        self.coi_blast_db = coi_blast_db
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

            Args:
                None

            Returns:
                String: The output to the taxonomy.sqlite database
        """
        PathManager.mkdir_p(os.path.join(self.coi_blast_db))

        ####################################################################
        #
        ####################################################################

        basename = 'map_taxids.tsv'
        map_taxids_tsv_url = os.path.join(public_data_dir, basename)
        map_taxids_tsv_path = os.path.join(self.coi_blast_db, basename)
        if not os.path.isfile(os.path.join(map_taxids_tsv_path)):
            urllib.request.urlretrieve(map_taxids_tsv_url, map_taxids_tsv_path)

        ####################################################################
        #
        ####################################################################

        basename = 'nt_coi.nhr'
        coi_db_url = os.path.join(public_data_dir, basename)
        coi_db_path = os.path.join(self.coi_blast_db, basename)
        if not os.path.isfile(os.path.join(coi_db_path)):
            urllib.request.urlretrieve(coi_db_url, coi_db_path)

        ####################################################################
        #
        ####################################################################

        basename = 'nt_coi.nin'
        coi_db_url = os.path.join(public_data_dir, basename)
        coi_db_path = os.path.join(self.coi_blast_db, basename)
        if not os.path.isfile(os.path.join(coi_db_path)):
            urllib.request.urlretrieve(coi_db_url, coi_db_path)

        ####################################################################
        #
        ####################################################################

        basename = 'nt_coi.nog'
        coi_db_url = os.path.join(public_data_dir, basename)
        coi_db_path = os.path.join(self.coi_blast_db, basename)
        if not os.path.isfile(os.path.join(coi_db_path)):
            urllib.request.urlretrieve(coi_db_url, coi_db_path)

        ####################################################################
        #
        ####################################################################

        basename = 'nt_coi.nsd'
        coi_db_url = os.path.join(public_data_dir, basename)
        coi_db_path = os.path.join(self.coi_blast_db, basename)
        if not os.path.isfile(os.path.join(coi_db_path)):
            urllib.request.urlretrieve(coi_db_url, coi_db_path)

        ####################################################################
        #
        ####################################################################

        basename = 'nt_coi.nsi'
        coi_db_url = os.path.join(public_data_dir, basename)
        coi_db_path = os.path.join(self.coi_blast_db, basename)
        if not os.path.isfile(os.path.join(coi_db_path)):
            urllib.request.urlretrieve(coi_db_url, coi_db_path)

        ####################################################################
        #
        ####################################################################

        basename = "nt_coi.nsq"
        coi_db_url = os.path.join(public_data_dir, basename)
        coi_db_path = os.path.join(self.coi_blast_db, basename)
        if not os.path.isfile(os.path.join(coi_db_path)):
            urllib.request.urlretrieve(coi_db_url, coi_db_path)

        ####################################################################
        #
        ####################################################################

    @staticmethod
    def create_parser():
        parser = argparse.ArgumentParser()
        parser.add_argument('--coi_blast_db', dest='coi_blast_db', action='store', help="Path COI Blast DB",
                            required=True)
        return parser

    @classmethod
    def main(cls):
        parser = DBblastCOI.create_parser()
        args = parser.parse_args()
        coi_blast_db = DBblastCOI(coi_blast_db=vars(args)['coi_blast_db'])
        coi_blast_db.download()
