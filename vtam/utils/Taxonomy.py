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

    def get_one_tax_id_lineage(self, tax_id):
        """
        Takes a tax_id and creates a dictionary with the taxonomy lineage in
        this form {'species': 183142, 'genus': 10194, 'family': 10193, 'order': 84394,
     'superorder': 1709201, 'class': 10191, 'phylum': 10190, 'no rank': 131567, 'kingdom': 33208,
     'superkingdom': 2759}

        Parameters
        ----------
        tax_id : int
            NCBI taxon id

        Returns
        -------
        dic
        Dictionnary with taxonomy lineage for given tax_id

        """

        # lineage_dic = {'tax_id': tax_id}
        lineage_dic = {}
        while tax_id != 1:
            # tax_id is found as normal tax in the taxonomy file
            if tax_id in self.df.index:
                tax_id_row = self.df.loc[tax_id, ]
            # tax_id is found as old_tax_id column in the taxonomy file
            elif tax_id in self.old_tax_df.index.tolist():  # Try old tax id
                tax_id_new = self.old_tax_df.loc[tax_id, 'tax_id']
                tax_id_row = self.df.loc[tax_id_new, ]
            # tax_id is not found in the taxonomy file.
            # Return current lineage dic and exit the function
            else:
                Logger.instance().warning(
                    "The taxon ID {} in the Blast database is missing in the taxonomy.tsv. "
                    "Consider updating this file with the following command: vtam taxonomy --output taxonomy.tsv.".format(tax_id))
                # raise VTAMexception("tax_id {} from Blast database not found in the taxonomy.tsv file".format(tax_id))
                return lineage_dic
            rank = tax_id_row['rank']
            parent_tax_id = tax_id_row['parent_tax_id']
            lineage_dic[rank] = tax_id
            tax_id = parent_tax_id
        return lineage_dic

    def get_several_tax_id_lineages(self, tax_id_list):
        """
        Takes a list of tax_id's and creates a DataFrame with the taxonomy lineages in columns
        and the tax_id as index

tax_id (index)  no rank    species     genus     family     order     class
1246992   131567   741276.0    5533.0  1799696.0  231213.0  162481.0    29000.0
1112827   131567  1112827.0    6220.0   941271.0    6219.0    6218.0

        Parameters
        ----------
        tax_id : int
            NCBI taxon id

        Returns
        -------
        DataFrame
        DataFrame with lineages in columns and tax_id as index

        """

        lineage_list = []
        for target_tax_id_i, target_tax_id in enumerate(tax_id_list):
            if target_tax_id_i % 100 == 0:
                Logger.instance().debug(
                    "Get lineage of {}-th tax id {} (Total {} tax ids)".format(
                        target_tax_id_i, target_tax_id, len(tax_id_list)))
            lineage_list.append({**{'tax_id': target_tax_id},
                                 **self.get_one_tax_id_lineage(tax_id=target_tax_id)})
        tax_id_lineage_df = pandas.DataFrame(lineage_list)
        tax_id_lineage_df.set_index('tax_id', drop=True, inplace=True, verify_integrity=True)
        return tax_id_lineage_df
