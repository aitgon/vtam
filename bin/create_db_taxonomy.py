#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import inspect
import urllib.request
import tarfile

import errno, pandas, urllib

import os

from wopmetabarcoding.utils.constants import data_dir
from wopmetabarcoding.utils.logger import logger

from sqlalchemy import create_engine

def download():
    try:
        os.makedirs(data_dir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    # Download files
    remotefile = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
    new_taxdump_path = os.path.join(data_dir, os.path.basename(remotefile))
    logger.debug(
        "file: {}; line: {}; Downloading NCBI taxonomy dump".format(__file__, inspect.currentframe().f_lineno))
    if not os.path.isfile(new_taxdump_path): # TODO verify MD5
        urllib.request.urlretrieve(remotefile, new_taxdump_path)
    return new_taxdump_path


def f_create_taxonomy_db(path_taxonomy_db_sqlite):
    new_taxdump_path = download()
    logger.debug(
        "file: {}; line: {}; Extracting NCBI taxonomy dump".format(__file__, inspect.currentframe().f_lineno))
    if not (os.path.isfile(os.path.join(os.path.dirname(new_taxdump_path), "nodes.dmp"))\
          and os.path.isfile(os.path.join(os.path.dirname(new_taxdump_path), "names.dmp"))\
          and os.path.isfile(os.path.join(os.path.dirname(new_taxdump_path), "merged.dmp"))): # TODO verify MD5
        tar = tarfile.open(new_taxdump_path, "r:gz")
        tar.extractall(path=data_dir)
        tar.close()
    logger.debug(
        "file: {}; line: {}; Reading and processing NCBI taxonomy dump".format(__file__, inspect.currentframe().f_lineno))
    #
    nodes_dmp = os.path.join(data_dir, "nodes.dmp")
    nodes_dmp_df = pandas.read_table(nodes_dmp, header=None, sep='\t', engine='python', usecols=[0, 2, 4],
                      names=['tax_id', 'parent_tax_id', 'rank'])
    #
    names_dmp = os.path.join(data_dir, "names.dmp")
    names_dmp_df = pandas.read_table(names_dmp, header=None, sep=r'\t', engine='python', usecols=[0, 2, 6],
                                     names=['tax_id', 'name_txt', 'name_class'])
    names_dmp_df = names_dmp_df.loc[names_dmp_df.name_class=='scientific name']
    names_dmp_df = names_dmp_df[['tax_id', 'name_txt']]
    #
    taxonomy_df = nodes_dmp_df.merge(names_dmp_df, on='tax_id')
    # del(nodes_dmp_df)
    # del(names_dmp_df)
    #
    merged_dmp = os.path.join(data_dir, "merged.dmp")
    merged_dmp_df = pandas.read_table(merged_dmp, header=None, sep='\t', engine='python', usecols=[0, 2],
                                      names=['old_tax_id', 'tax_id'])
    #
    taxonomy_df = taxonomy_df.merge(merged_dmp_df, on='tax_id', how='left')
    #
    logger.debug(
        "file: {}; line: {}; Write to sqlite DB".format(__file__, inspect.currentframe().f_lineno))
    engine = create_engine('sqlite:///{}'.format(path_taxonomy_db_sqlite), echo=False)
    try:
        taxonomy_df.to_sql('taxonomy', con=engine, index = False)
    except ValueError:
        return path_taxonomy_db_sqlite
    #
    return path_taxonomy_db_sqlite

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', dest='output', help="Path to sqlite DB with taxonomy information")
    return parser

def main():
    parser = create_parser()
    args = parser.parse_args()
    f_create_taxonomy_db(args.output)

if __name__=='__main__':
    main()
