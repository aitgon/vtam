from vtam.utils.VTAMexception import VTAMexception
from vtam.utils.Logger import Logger

import inspect
import os
import pandas
import sqlalchemy
import tarfile
import urllib
import urllib.request

from vtam.utils.PathManager import PathManager
from vtam.utils.constants import url_taxonomy_tsv
from sqlalchemy.exc import OperationalError

class DBtaxonomy(object):

    def __init__(self, output=None, precomputed=True):
        #
        # Download the precomputed database. The alternative to create a new DB with the create_db_taxonomy executable
        self.precomputed = precomputed
        #
        # output to the taxonomy.sqlite file
        self.output = os.path.join(os.getcwd(), "taxonomy.tsv")
        if not output is None:
            self.output = output
        #
        self.tempdir = PathManager.instance().get_tempdir()

    def get_path(self):
        if not os.path.isfile(self.output):
            if self.precomputed:
                self.download_taxonomy_tsv()
            else:
                self.create_taxonomy_db()
        return self.output


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
              and os.path.isfile(os.path.join(os.path.dirname(new_taxdump_path), "merged.dmp"))):
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
        # engine = sqlalchemy.create_engine('sqlite:///{}'.format(self.output), echo=False)
        try:
            taxonomy_df.to_csv(self.output, sep="\t", header=True, float_format='%.0f', index=False)
            # taxonomy_df.to_sql('taxonomy', con=engine, index = False)
        except ValueError as valerr:
            Logger.instance().error(VTAMexception("{}. Error during the creation of the taxonomy DB".format(valerr)))
        except sqlalchemy.exc.OperationalError as opererr:
            Logger.instance().error(
                VTAMexception("{}. Please, verify the output argument: {}".format(opererr, self.output)))

    ##########################################################
    #
    # Define/create taxonomy.tsv output
    #
    ##########################################################
    def download_taxonomy_tsv(self):
        """
        Copy the online SQLITE taxonomy DB at "http://pedagogix-tagc.univ-mrs.fr/~gonzalez/vtam/taxonomy.sqlite"
        to the pathname output
        """
        Logger.instance().debug(
            "file: {}; line: {}; download_taxonomy_tsv()".format(__file__,
                                                                    inspect.currentframe().f_lineno, ))
        if not os.path.isfile(self.output):
            urllib.request.urlretrieve(url_taxonomy_tsv, self.output)

    # @staticmethod
    # def create_parser():
    #     parser = argparse.ArgumentParser()
    #     parser.add_argument('-o', '--output', dest='output', action='store', help="Path to custom COI blast db",
    #                         required=True)
    #     parser.add_argument('--precomputed', dest='precomputed', action='store_true', default=True,
    #                         help="Path to custom COI blast db",
    #                         required=False)
    #     return parser
    #
    # @classmethod
    # def main(cls):
    #     parser = DBtaxonomy.create_parser()
    #     args = parser.parse_args()
    #     taxonomydb = DBtaxonomy(output=vars(args)['output'], precomputed=vars(args)['precomputed'], )
    #     taxonomydb.create_taxonomy_db()

