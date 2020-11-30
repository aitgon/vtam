import inspect
import os
import pandas
import pathlib

from vtam.utils.FileParams import FileParams

from vtam.utils.Logger import Logger
from vtam.utils.PathManager import PathManager
from vtam.utils.RunnerBlast import RunnerBlast
from vtam.utils.RunnerLTGselection import RunnerLTGselection
from vtam.utils.VTAMexception import VTAMexception
from vtam.utils.constants import rank_hierarchy, identity_list


class RunnerTaxAssign(object):
    """Will assign variants to a taxon"""

    def __init__(self, sequence_list, taxonomy_df, blast_db_dir, blast_db_name,
             num_threads, params):

        # self.variant_df = variant_df
        # stores tax_id and old_tax_id
        self.old_tax_id_df = taxonomy_df[['old_tax_id']].drop_duplicates()
        self.taxonomy_df = taxonomy_df[[
            'parent_tax_id', 'rank', 'name_txt']].drop_duplicates()
        self.blast_db_dir = blast_db_dir
        self.this_temp_dir = os.path.join(PathManager.instance().get_tempdir(),
            os.path.basename(__file__))
        pathlib.Path(self.this_temp_dir).mkdir(exist_ok=True)

        self.num_threads = num_threads

        #######################################################################
        #
        # Parameters
        #
        #######################################################################

        params_dic = FileParams(params).get_params_dic()
        qcov_hsp_perc = params_dic['qcov_hsp_perc']

        self.ltg_rule_threshold = params_dic['ltg_rule_threshold']
        self.include_prop = params_dic['include_prop']
        self.min_number_of_taxa = params_dic['min_number_of_taxa']

        #######################################################################
        #
        # 2 Create FASTA file with Variants
        #
        #######################################################################

        Logger.instance().debug(
            "file: {}; line: {}; Create SortedReadFile from Variants".format(
                __file__, inspect.currentframe().f_lineno))
        variant_fasta = os.path.join(self.this_temp_dir, 'variant.fasta')
        # variant_df_utils = DataframeVariant(variant_df)
        # variant_df_utils.to_fasta(variant_fasta)
        with open(variant_fasta, 'w') as fout:
            for seq in sequence_list:
                fout.write(">{}\n{}\n".format(seq, seq))

        #######################################################################
        #
        # 3 Run local blast
        #
        #######################################################################

        runner_blast = RunnerBlast(variant_fasta, blast_db_dir, blast_db_name,
            num_threads, qcov_hsp_perc)
        # run blast
        blast_output_tsv = runner_blast.run_local_blast()
        # process blast results
        blast_output_df = RunnerBlast.process_blast_result(blast_output_tsv)

        self.ltg_df = None  # Init it

        #######################################################################
        #
        # Read target_tax_id
        # Compute lineages for each unique target_tax_id
        # Create a DF with these columns: tax_id and its lineage in wide format
        #
        #######################################################################

        Logger.instance().debug("file: {}; line: {}; Open taxonomy.tsv DB".format(
                __file__, inspect.currentframe().f_lineno))
        blast_output_df.target_tax_id = pandas.to_numeric(blast_output_df.target_tax_id)
        #
        Logger.instance().debug(
            "file: {}; line: {}; Annotate each target_tax_id with its lineage as columns in wide format".format(
                __file__, inspect.currentframe().f_lineno))
        lineage_list = []
        for target_tax_id_i, target_tax_id in enumerate(
                blast_output_df.target_tax_id.unique().tolist()):
            if target_tax_id_i % 100 == 0:
                Logger.instance().debug(
                    "Get lineage of {}-th tax id {} (Total {} tax ids)".format(
                        target_tax_id_i, target_tax_id, len(
                            blast_output_df.target_tax_id.unique().tolist())))
            lineage_list.append(
                self.tax_id_to_taxonomy_lineage(
                    tax_id=target_tax_id))
        tax_id_to_lineage_df = pandas.DataFrame(lineage_list)

        #######################################################################
        #
        # Merge tax lineages and the blast result
        #
        #######################################################################

        Logger.instance().debug(
            "file: {}; line: {}; Merge blast result including tax_id with their lineages".format(
                __file__, inspect.currentframe().f_lineno))
        # Merge local blast output with tax_id_to_lineage_df
        variantid_identity_lineage_df = blast_output_df.merge(
            tax_id_to_lineage_df, left_on='target_tax_id', right_on='tax_id')
        variantid_identity_lineage_df.drop('tax_id', axis=1, inplace=True)

        """(Pdb) variantid_identity_lineage_df.columns
Index(['variant_id', 'target_id', 'identity', 'evalue', 'coverage',
       'target_tax_id', 'no rank', 'species', 'genus', 'family', 'order',
       'class', 'subphylum', 'phylum', 'subkingdom', 'kingdom', 'superkingdom',
       'superfamily', 'infraorder', 'suborder', 'infraclass', 'subclass',
       'tribe', 'subfamily', 'cohort', 'subgenus', 'subspecies', 'parvorder',
       'superorder', 'subcohort', 'superclass', 'species group', 'subtribe',
       'section', 'varietas', 'species subgroup'],
      dtype='object')"""

        #######################################################################
        #
        #  blast_output_to_ltg_tax_id
        # this function returns a data frame containing the Ltg rank and Ltg Tax_id for each variant
        #
        #######################################################################

        Logger.instance().debug(
            "file: {}; line: {}; Main loop over variant and identity to"
            "compute the whole set of ltg_tax_id and ltg_rank for each variant_id"
            "to a dataframe".format(
                __file__, inspect.currentframe().f_lineno))
        # self.ltg_df = self.blast_output_to_ltg_tax_id(variantid_identity_lineage_df)
        runner_ltg_selection = RunnerLTGselection(variantid_identity_lineage_df=variantid_identity_lineage_df,
                                                  taxonomy_df=self.taxonomy_df, params=params)
        self.ltg_df = runner_ltg_selection.blast_output_to_ltg_tax_id()
        # import pdb; pdb.set_trace()

    def tax_id_to_taxonomy_lineage(self, tax_id):
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

        # # check if tax_id in taxonomy_df or old_tax_id and if not raise error and exit
        # if not (tax_id in self.taxonomy_df.index or tax_id in self.old_tax_id_df.old_tax_id.tolist()):
        #     Logger.instance().error(
        #         "The taxon ID {} in the Blast database is missing in the taxonomy.tsv. "
        #         "Consider updating this file with the following command: vtam taxonomy --output taxonomy.tsv. "
        #         "The workflow will exit.".format(tax_id))
        #     raise VTAMexception("tax_id {} from Blast database not found in the taxonomy.tsv file".format(tax_id))

        lineage_dic = {'tax_id': tax_id}
        while tax_id != 1:
            # tax_id is found as normal index in the taxonomy file
            if tax_id in self.taxonomy_df.index:
                tax_id_row = self.taxonomy_df.loc[tax_id, ]
            # tax_id is found as old_tax_id column in the taxonomy file
            elif tax_id in self.old_tax_id_df.old_tax_id.tolist():  # Try old tax id
                tax_id2 = self.old_tax_id_df.loc[self.old_tax_id_df.old_tax_id ==
                                                 tax_id, 'old_tax_id'].index[0]
                tax_id_row = self.taxonomy_df.loc[tax_id2, ]
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
