from pathlib import Path

from vtam import VTAMexception
from vtam.utils.Logger import Logger

import inspect
import os
import pandas
import sqlalchemy
import tarfile
import urllib
import urllib.request

from vtam.utils.PathManager import PathManager
from vtam.utils.constants import url_taxonomy_sqlite


class TaxonomyDB(object):

    def __init__(self, path=None, package=True):
        #
        # Download the package database. The alternative to create a new DB with the create_db_taxonomy executable
        self.package = package
        #
        # path to the taxonomy.sqlite file
        self.path = os.path.join(os.getcwd(), "taxonomy.sqlite")
        if not path is None:
            self.path = path
        #
        self.tempdir = PathManager.instance().get_tempdir()

    def get_path(self):
        if not os.path.isfile(self.path):
            if self.package:
                self.__download_taxonomy_sqlite()
            else:
                self.create_taxonomy_db()
        return self.path


    def __download_ncbi_taxonomy_dump(self):
        # Download files
        remotefile = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
        new_taxdump_path = os.path.join(self.tempdir, os.path.basename(remotefile))
        Logger.instance().debug(
            "file: {}; line: {}; Downloading NCBI taxonomy dump".format(__file__, inspect.currentframe().f_lineno))
        if not os.path.isfile(new_taxdump_path):
            urllib.request.urlretrieve(remotefile, new_taxdump_path)
        return new_taxdump_path

    def create_taxonomy_db(self):
        new_taxdump_path = self.__download_ncbi_taxonomy_dump()
        #
        Logger.instance().debug(
            "file: {}; line: {}; Extracting NCBI taxonomy dump".format(__file__, inspect.currentframe().f_lineno))
        if not (os.path.isfile(os.path.join(os.path.dirname(new_taxdump_path), "nodes.dmp"))\
              and os.path.isfile(os.path.join(os.path.dirname(new_taxdump_path), "names.dmp"))\
              and os.path.isfile(os.path.join(os.path.dirname(new_taxdump_path), "merged.dmp"))): # TODO verify MD5
            tar = tarfile.open(new_taxdump_path, "r:gz")
            tar.extractall(path=self.tempdir)
            tar.close()
        Logger.instance().debug(
            "file: {}; line: {}; Reading and processing NCBI taxonomy dump".format(__file__, inspect.currentframe().f_lineno))
        #
        nodes_dmp = os.path.join(self.tempdir, "nodes.dmp")
        nodes_dmp_df = pandas.read_table(nodes_dmp, header=None, sep='\t', engine='python', usecols=[0, 2, 4],
                          names=['tax_id', 'parent_tax_id', 'rank'])
        #
        names_dmp = os.path.join(self.tempdir, "names.dmp")
        names_dmp_df = pandas.read_table(names_dmp, header=None, sep=r'\t', engine='python', usecols=[0, 2, 6],
                                         names=['tax_id', 'name_txt', 'name_class'])
        names_dmp_df = names_dmp_df.loc[names_dmp_df.name_class=='scientific name']
        names_dmp_df = names_dmp_df[['tax_id', 'name_txt']]
        #
        taxonomy_df = nodes_dmp_df.merge(names_dmp_df, on='tax_id')
        #
        merged_dmp = os.path.join(self.tempdir, "merged.dmp")
        merged_dmp_df = pandas.read_table(merged_dmp, header=None, sep='\t', engine='python', usecols=[0, 2],
                                          names=['old_tax_id', 'tax_id'])
        #
        taxonomy_df = taxonomy_df.merge(merged_dmp_df, on='tax_id', how='left')
        #
        Logger.instance().debug(
            "file: {}; line: {}; Write to sqlite DB".format(__file__, inspect.currentframe().f_lineno))
        engine = sqlalchemy.create_engine('sqlite:///{}'.format(self.path), echo=False)
        try:
            taxonomy_df.to_sql('taxonomy', con=engine, index = False)
        except ValueError as valerr:
            Logger.instance().error(VTAMexception("Error during the creation of the taxonomy DB"))
        # #
        # return path_taxonomy_db_sqlite_path

    ##########################################################
    #
    # Define/create taxonomy.sqlite path
    #
    ##########################################################
    def __download_taxonomy_sqlite(self):
        """
        Copy the online SQLITE taxonomy DB at "http://pedagogix-tagc.univ-mrs.fr/~gonzalez/vtam/taxonomy.sqlite"
        to the pathname path
        """
        Logger.instance().debug(
            "file: {}; line: {}; __download_taxonomy_sqlite()".format(__file__,
                                                                    inspect.currentframe().f_lineno, ))
        if not os.path.isfile(self.path):
            urllib.request.urlretrieve(url_taxonomy_sqlite, self.path)
