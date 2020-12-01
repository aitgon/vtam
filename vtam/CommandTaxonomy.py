import gzip
import inspect
import os
import pandas
import pathlib
import shutil
import sqlalchemy
import tarfile
import urllib
import urllib.request

from sqlalchemy.exc import OperationalError

from vtam.utils import tqdm_hook
from vtam.utils.Logger import Logger
from vtam.utils.PathManager import PathManager
from vtam.utils.VTAMexception import VTAMexception
from vtam.utils.constants import taxonomy_tsv_gz_url1, taxonomy_tsv_gz_url2, taxonomy_tsv_gz_url3
from tqdm import tqdm


class CommandTaxonomy(object):

    def __init__(self, taxonomy_tsv=None):
        """

        :param taxonomy_tsv: Path to the taxonomy_tsv. Default None
        :type taxonomy_tsv: str

        :rtype: None
        """

        if taxonomy_tsv is None:  # If None, download to current wdir
            self.taxonomy_tsv_path = os.path.join(os.getcwd(), "taxonomy.tsv")
        else:  # Download to tsv_path
            self.taxonomy_tsv_path = taxonomy_tsv

        pathlib.Path(
            os.path.dirname(taxonomy_tsv)).mkdir(
            parents=True, exist_ok=True)

        self.tempdir = PathManager.instance().get_tempdir()

        package_path = os.path.join(PathManager.get_package_path())
        self.taxonomy_tsv_gz_path = os.path.join(package_path, "..", "data", "taxonomy.tsv.gz")

    def __download_ncbi_taxonomy_dump(self):
        # Download files
        remotefile = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
        new_taxdump_path = os.path.join(
            self.tempdir, os.path.basename(remotefile))
        Logger.instance().debug(
            "file: {}; line: {}; Downloading NCBI taxonomy dump".format(
                __file__, inspect.currentframe().f_lineno))
        if not os.path.isfile(new_taxdump_path):
            Logger.instance().info("Downloading NCBI taxonomy dump")
            # urllib.request.urlretrieve(remotefile, new_taxdump_path, MyProgressBar())
            with tqdm(...) as t:
                t.set_description(os.path.basename(new_taxdump_path))
                urllib.request.urlretrieve(remotefile, new_taxdump_path, reporthook=tqdm_hook(t))
        return new_taxdump_path

    def create_denovo_from_ncbi(self):
        new_taxdump_path = self.__download_ncbi_taxonomy_dump()
        #
        Logger.instance().debug(
            "file: {}; line: {}; Extracting NCBI taxonomy dump".format(
                __file__, inspect.currentframe().f_lineno))
        if not (
            os.path.isfile(
                os.path.join(
                    os.path.dirname(new_taxdump_path),
                    "nodes.dmp")) and os.path.isfile(
                os.path.join(
                    os.path.dirname(new_taxdump_path),
                    "names.dmp")) and os.path.isfile(
                        os.path.join(
                            os.path.dirname(new_taxdump_path),
                            "merged.dmp"))):
            tar = tarfile.open(new_taxdump_path, "r:gz")
            tar.extractall(path=self.tempdir)
            tar.close()
        Logger.instance().debug(
            "file: {}; line: {}; Reading and processing NCBI taxonomy dump".format(
                __file__, inspect.currentframe().f_lineno))
        #
        nodes_dmp = os.path.join(self.tempdir, "nodes.dmp")
        nodes_dmp_df = pandas.read_table(
            nodes_dmp, header=None, sep='\t', engine='python', usecols=[
                0, 2, 4], names=[
                'tax_id', 'parent_tax_id', 'rank'])
        #
        names_dmp = os.path.join(self.tempdir, "names.dmp")
        names_dmp_df = pandas.read_table(
            names_dmp, header=None, sep=r'\t', engine='python', usecols=[
                0, 2, 6], names=[
                'tax_id', 'name_txt', 'name_class'])
        names_dmp_df = names_dmp_df.loc[names_dmp_df.name_class ==
                                        'scientific name']
        names_dmp_df = names_dmp_df[['tax_id', 'name_txt']]
        #
        taxonomy_df = nodes_dmp_df.merge(names_dmp_df, on='tax_id')
        #
        merged_dmp = os.path.join(self.tempdir, "merged.dmp")
        merged_dmp_df = pandas.read_table(
            merged_dmp, header=None, sep='\t', engine='python', usecols=[
                0, 2], names=[
                'old_tax_id', 'tax_id'])
        #
        taxonomy_df = taxonomy_df.merge(merged_dmp_df, on='tax_id', how='left')
        #
        Logger.instance().debug(
            "file: {}; line: {}; Write to TSV DB".format(
                __file__, inspect.currentframe().f_lineno))
        try:
            taxonomy_df.to_csv(
                self.taxonomy_tsv_path,
                sep="\t",
                header=True,
                float_format='%.0f',
                index=False)
        except ValueError as valerr:
            Logger.instance().error(
                VTAMexception(
                    "{}. Error during the creation of the taxonomy DB".format(valerr)))
        except sqlalchemy.exc.OperationalError as opererr:
            Logger.instance().error(
                VTAMexception(
                    "{}. Please, verify the output argument: {}".format(
                        opererr, self.taxonomy_tsv_path)))

    def download_precomputed_taxonomy(self):
        """
        Copy the online TSV taxonomy DB
        to the pathname output
        """
        Logger.instance().debug(
            "file: {}; line: {}; Downloading taxonomy tsv".format(
                __file__, inspect.currentframe().f_lineno, ))

        ############################################################################################
        #
        # Download sorted reads dataset
        #
        ############################################################################################

        taxonomy_tsv_gz_path = '{}.gz'.format(self.taxonomy_tsv_path)
        # Test first in local dir, otherwise in the remote URLs
        if not os.path.isfile(self.taxonomy_tsv_path) or pathlib.Path(self.taxonomy_tsv_path).stat().st_size < 1000000:
            try:
                # urllib.request.urlretrieve(taxonomy_tsv_gz_url1, taxonomy_tsv_gz_path, MyProgressBar())
                with tqdm(...) as t:
                    t.set_description(os.path.basename(taxonomy_tsv_gz_url1))
                    urllib.request.urlretrieve(taxonomy_tsv_gz_url1, taxonomy_tsv_gz_path, reporthook=tqdm_hook(t))
            except Exception:
                try:
                    # urllib.request.urlretrieve(taxonomy_tsv_gz_url2, taxonomy_tsv_gz_path,
                    #                            MyProgressBar())
                    with tqdm(...) as t:
                        t.set_description(os.path.basename(taxonomy_tsv_gz_url2))
                        urllib.request.urlretrieve(taxonomy_tsv_gz_url2, taxonomy_tsv_gz_path, reporthook=tqdm_hook(t))
                except Exception:
                    # urllib.request.urlretrieve(taxonomy_tsv_gz_url3, taxonomy_tsv_gz_path,
                    #                            MyProgressBar())
                    with tqdm(...) as t:
                        t.set_description(os.path.basename(taxonomy_tsv_gz_url1))
                        urllib.request.urlretrieve(taxonomy_tsv_gz_url3, taxonomy_tsv_gz_path, reporthook=tqdm_hook(t))
            with gzip.open(taxonomy_tsv_gz_path, 'rb') as fin:
                with open(self.taxonomy_tsv_path, 'wb') as fout:
                    shutil.copyfileobj(fin, fout)
            try:
                pathlib.Path(taxonomy_tsv_gz_path).unlink()
            except FileNotFoundError:
                pass

    def main(self, precomputed=True):
        """

        :param precomputed: Download precomputed taxonomy file. Default False
        :type precomputed: bool

        :rtype: None
        """

        if precomputed:
            self.download_precomputed_taxonomy()
        else:
            self.create_denovo_from_ncbi()
