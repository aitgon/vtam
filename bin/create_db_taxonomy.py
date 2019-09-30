#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import inspect
import tempfile
import urllib.request
import tarfile

import errno, pandas, urllib

import os

from wopmetabarcoding.utils.Logger import Logger

from sqlalchemy import create_engine


class TaxonomyDBCreator(object):


    def __init__(self, datadir=None):
        self.datadir = datadir
        if self.datadir is None:
            self.datadir = tempfile.mkdtemp()

    def download(self):
        # Download files
        remotefile = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
        new_taxdump_path = os.path.join(self.datadir, os.path.basename(remotefile))
        Logger.instance().debug(
            "file: {}; line: {}; Downloading NCBI taxonomy dump".format(__file__, inspect.currentframe().f_lineno))
        if not os.path.isfile(new_taxdump_path):
            urllib.request.urlretrieve(remotefile, new_taxdump_path)
        return new_taxdump_path


    def f_create_taxonomy_db(self, path_taxonomy_db_sqlite_path, datadir=None):
        new_taxdump_path = self.download()
        Logger.instance().debug(
            "file: {}; line: {}; Extracting NCBI taxonomy dump".format(__file__, inspect.currentframe().f_lineno))
        if not (os.path.isfile(os.path.join(os.path.dirname(new_taxdump_path), "nodes.dmp"))\
              and os.path.isfile(os.path.join(os.path.dirname(new_taxdump_path), "names.dmp"))\
              and os.path.isfile(os.path.join(os.path.dirname(new_taxdump_path), "merged.dmp"))): # TODO verify MD5
            tar = tarfile.open(new_taxdump_path, "r:gz")
            tar.extractall(path=self.datadir)
            tar.close()
        Logger.instance().debug(
            "file: {}; line: {}; Reading and processing NCBI taxonomy dump".format(__file__, inspect.currentframe().f_lineno))
        #
        nodes_dmp = os.path.join(self.datadir, "nodes.dmp")
        nodes_dmp_df = pandas.read_table(nodes_dmp, header=None, sep='\t', engine='python', usecols=[0, 2, 4],
                          names=['tax_id', 'parent_tax_id', 'rank'])
        #
        names_dmp = os.path.join(self.datadir, "names.dmp")
        names_dmp_df = pandas.read_table(names_dmp, header=None, sep=r'\t', engine='python', usecols=[0, 2, 6],
                                         names=['tax_id', 'name_txt', 'name_class'])
        names_dmp_df = names_dmp_df.loc[names_dmp_df.name_class=='scientific name']
        names_dmp_df = names_dmp_df[['tax_id', 'name_txt']]
        #
        taxonomy_df = nodes_dmp_df.merge(names_dmp_df, on='tax_id')
        #
        merged_dmp = os.path.join(self.datadir, "merged.dmp")
        merged_dmp_df = pandas.read_table(merged_dmp, header=None, sep='\t', engine='python', usecols=[0, 2],
                                          names=['old_tax_id', 'tax_id'])
        #
        taxonomy_df = taxonomy_df.merge(merged_dmp_df, on='tax_id', how='left')
        #
        Logger.instance().debug(
            "file: {}; line: {}; Write to sqlite DB".format(__file__, inspect.currentframe().f_lineno))
        engine = create_engine('sqlite:///{}'.format(path_taxonomy_db_sqlite_path), echo=False)
        try:
            taxonomy_df.to_sql('taxonomy', con=engine, index = False)
        except ValueError:
            return path_taxonomy_db_sqlite_path
        #
        return path_taxonomy_db_sqlite_path


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', dest='output', action='store', help="Path to sqlite DB with taxonomy information",
                        required=True)
    parser.add_argument('--datadir', action='store',
            help="Directory that will be used for download NCBI data. If this argument is absent, the data will be downloaded to the temporary directory."
                        , required=False)
    return parser

def main():
    parser = create_parser()
    args = parser.parse_args()
    taxonomydbcreator = TaxonomyDBCreator(args.datadir)
    taxonomydbcreator.f_create_taxonomy_db(args.output)


if __name__=='__main__':
    main()
