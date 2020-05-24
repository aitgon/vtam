import pandas

from vtam.utils.constants import rank_hierarchy_asv_table


class TaxLineage(object):
    """This class construct a TaxLineage for a given tax_id and based on the taxonomic_tsv file"""

    def __init__(self, taxonomic_tsv_path):
        self.taxonomy_df = pandas.read_csv(
            taxonomic_tsv_path, sep="\t", header=0, dtype={
                'tax_id': 'int', 'parent_tax_id': 'int', 'old_tax_id': 'float'})

    def create_lineage_from_one_tax_id(self, tax_id, tax_name=False):
        """
        Takes tax_id and creates a dictionary with the taxonomy lineage

        Args:
            tax_id (Integer): Identifier of taxon
            tax_name (String): Store tax_name instead of tax_id

        Returns:
            Dictionnary: with taxonomy information for given tax_id,
            {'tax_id': 183142, 'species': 183142, 'genus': 10194, 'family': 10193, 'order': 84394,
                             'superorder': 1709201, 'class': 10191, 'phylum': 10190, 'no rank': 131567, 'kingdom': 33208,
                             'superkingdom': 2759}

        """

        # Try to convert to int or return None otherwise
        try:
            tax_id = int(tax_id)
        except ValueError:
            return None
        except TypeError:
            return None

        tax_lineage_dic = {}
        tax_lineage_dic['tax_id'] = tax_id
        while tax_id != 1:
            # try to use taxonomy_df.tax_id
            tax_id_row = self.taxonomy_df.loc[self.taxonomy_df.tax_id == tax_id, ]
            # if row empty, try to use old_tax_id
            if tax_id_row.shape[0] == 0:
                tax_id_row = self.taxonomy_df.loc[self.taxonomy_df.old_tax_id == tax_id, ]
                # if row still empty, return None
                if tax_id_row.shape[0] == 0:
                    return None

            rank = tax_id_row['rank'].values[0]
            parent_tax_id = tax_id_row['parent_tax_id'].values[0]
            tax_lineage_dic[rank] = tax_id
            if tax_name:  # return tax_name instead of tax_id
                tax_name = tax_id_row['name_txt'].values[0]
                tax_lineage_dic[rank] = tax_name
            tax_id = parent_tax_id

        return tax_lineage_dic

    def create_lineage_from_tax_id_list(self, tax_id_list, tax_name=False):
        """
        Takes tax_id and creates a dictionary with the taxonomy lineage

        Args:
            tax_id_list (List): List of tax_ids
            tax_name (String): Append tax_name to dictionary

        Returns:
            Dictionnary: with taxonomy information for given tax_id,
            {'tax_id': 183142, 'species': 183142, 'genus': 10194, 'family': 10193, 'order': 84394,
                             'superorder': 1709201, 'class': 10191, 'phylum': 10190, 'no rank': 131567, 'kingdom': 33208,
                             'superkingdom': 2759}

        """
        taxa_lineage_list = []

        for tax_id in tax_id_list:

            # if not math.isnan(tax_id):
            # tax_lineage = TaxLineage(taxonomic_tsv_path=taxonomy_tsv_path)
            tax_lineage_dic = self.create_lineage_from_one_tax_id(
                tax_id=tax_id, tax_name=tax_name)
            # tax_lineage_dic = create_lineage_from_one_tax_id(tax_id, taxonomy_df, tax_name=True)
            if not (tax_lineage_dic is None):
                taxa_lineage_list.append(tax_lineage_dic)

        tax_lineage_df = pandas.DataFrame(data=taxa_lineage_list)
        # lineage_list_df_columns_sorted = list(
        #     filter(lambda x: x in lineage_df.columns.tolist(), rank_hierarchy_asv_table))
        # lineage_list_df_columns_sorted = lineage_list_df_columns_sorted + ['tax_id']
        # lineage_df = lineage_df[lineage_list_df_columns_sorted]

        lineage_list_df_columns_sorted = list(
            filter(
                lambda x: x in tax_lineage_df.columns.tolist(),
                rank_hierarchy_asv_table))
        lineage_list_df_columns_sorted = lineage_list_df_columns_sorted + \
            ['tax_id']
        tax_lineage_df = tax_lineage_df[lineage_list_df_columns_sorted]
        # do not move. required because sometimes tax_id is none
        tax_lineage_df = tax_lineage_df.astype({'tax_id': 'object'})

        return tax_lineage_df
