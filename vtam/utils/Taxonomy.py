import pandas
from vtam.utils.VTAMexception import VTAMexception

from vtam.utils.Logger import Logger


class Taxonomy(object):
    """A class for the taxonomy file"""

    def __init__(self, tsv=None, df=None):
        """Taxonomy gets initialize either from a TSV path or a DataFrame"""

        self.df = df
        if not (tsv is None):
            self.df = pandas.read_csv(tsv, sep="\t", header=0, dtype={'tax_id': 'int', 'parent_tax_id': 'int', 'old_tax_id': 'float'}).drop_duplicates()

        self.old_tax_df = self.df[['tax_id', 'old_tax_id']].drop_duplicates()
        self.old_tax_df = self.old_tax_df.loc[~self.old_tax_df.old_tax_id.isna()]
        self.old_tax_df.old_tax_id = self.old_tax_df.old_tax_id.astype('int')
        self.old_tax_df.set_index('old_tax_id', drop=True, inplace=True, verify_integrity=False)

        self.df = self.df.drop(['old_tax_id'], axis=1, inplace=False).drop_duplicates()
        self.df.set_index('tax_id', drop=True, inplace=True, verify_integrity=True)



    def create_lineage(self, tax_id):
        """
        Takes a tax_id and creates a dictionary with the taxonomy lineage

        :param tax_id: Identifier of taxon
        :type tax_id: int
        :param return_tax_name: (pandas.DataFrame: DataFrame with taxonomy information
        :type return_tax_name: str

        :return: Dic with ranks/keys and tax_id/values
        :rtype: list
        Args:
            tax_id (Integer): Identifier of taxon
            taxonomy_df (pandas.DataFrame: DataFrame with taxonomy information

        Returns:
            Dictionnary: with taxonomy information for given tax_id, {'tax_id': 183142, 'species': 183142, 'genus': 10194, 'family': 10193, 'order': 84394,
                             'superorder': 1709201, 'class': 10191, 'phylum': 10190, 'no rank': 131567, 'kingdom': 33208,
                             'superkingdom': 2759}

        """

        lineage_dic = {'tax_id': tax_id}
        while tax_id != 1:
            # tax_id is found as normal index in the taxonomy file
            if tax_id in self.df.index:
                tax_id_row = self.df.loc[tax_id, ]
            # tax_id is found as old_tax_id column in the taxonomy file
            elif tax_id in self.old_tax_df.old_tax_id.tolist():  # Try old tax id
                # tax_id2 = self.old_tax_id_df.loc[self.old_tax_id_df.old_tax_id ==
                #                                  tax_id, 'old_tax_id'].index[0]
                tax_id_new = self.old_tax_df.loc[tax_id, 'tax_id']
                tax_id_row = self.df.loc[tax_id_new, ]
            else:  # tax_id not in taxonomy, log error and exit
                Logger.instance().error(
                    "The taxon ID {} in the Blast database is missing in the taxonomy.tsv. "
                    "Consider updating this file with the following command: vtam taxonomy --output taxonomy.tsv. "
                    "The workflow will exit.".format(tax_id))
                raise VTAMexception("tax_id {} from Blast database not found in the taxonomy.tsv file".format(tax_id))
            rank = tax_id_row['rank']
            parent_tax_id = tax_id_row['parent_tax_id']
            lineage_dic[rank] = tax_id
            tax_id = parent_tax_id
        return lineage_dic
