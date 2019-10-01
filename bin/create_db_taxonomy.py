#!/usr/bin/env python
# -*- coding: utf-8 -*-

from wopmetabarcoding.utils.TaxonomyDB import TaxonomyDB

import argparse


import inspect
import tempfile
import urllib.request
import tarfile

import errno, pandas, urllib

import os

from wopmetabarcoding.utils.Logger import Logger

from sqlalchemy import create_engine





def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', dest='output', action='store', help="Path to sqlite DB with taxonomy information",
                        required=True)
    parser.add_argument('--datadir', action='store',
            help="Directory that will be used for __download_ncbi_taxonomy_dump NCBI data. If this argument is absent, the data will be downloaded to the temporary directory."
                        , required=False)
    return parser

def main():
    parser = create_parser()
    args = parser.parse_args()
    taxonomydb = TaxonomyDB(args.datadir)
    taxonomydb.create_taxonomy_db(args.output)


if __name__=='__main__':
    main()
