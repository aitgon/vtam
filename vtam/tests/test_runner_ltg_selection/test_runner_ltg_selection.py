import os
import pandas
import pathlib
import unittest

from vtam.utils.RunnerLTGselection import RunnerLTGselection
from vtam.utils.PathManager import PathManager
from vtam.CommandTaxonomy import CommandTaxonomy
from vtam.utils.Taxonomy import Taxonomy


class TestRunnerLTGselection(unittest.TestCase):

    def setUp(self):

        test_path = os.path.join(PathManager.get_test_path())
        self.variantid_identity_lineage_df = pandas.read_csv(os.path.join(test_path, "test_runner_ltg_selection/variantid_identity_lineage.tsv"), sep="\t", header=0)
        self.ltg_bak_df = pandas.read_csv(os.path.join(test_path, "test_runner_ltg_selection/ltg_bak.tsv"), sep="\t")

        # create_vtam_data_dir()
        testdir_path = os.path.join(PathManager.get_test_path())
        self.outdir_path = os.path.join(testdir_path, "outdir")
        pathlib.Path(self.outdir_path).mkdir(exist_ok=True, parents=True)
        taxonomy_tsv_path = os.path.join(self.outdir_path, "taxonomy.tsv")
        CommandTaxonomy(
            taxonomy_tsv=taxonomy_tsv_path).download_precomputed_taxonomy()

        self.taxonomy_df = pandas.read_csv(taxonomy_tsv_path, sep="\t", header=0,
                                      dtype={'tax_id': 'int', 'parent_tax_id': 'int',
                                             'old_tax_id': 'float'}).drop_duplicates()
        self.taxonomy_df.set_index('tax_id', drop=True, inplace=True)
        self.taxonomy_df = self.taxonomy_df[[
            'parent_tax_id', 'rank', 'name_txt']].drop_duplicates()
        taxonomy = Taxonomy(taxonomy_tsv_path)
        self.taxonomy_df = taxonomy.df

    def test_01(self):

        runner_ltg_selection = RunnerLTGselection(
            variantid_identity_lineage_df=self.variantid_identity_lineage_df, taxonomy_df=self.taxonomy_df, params=None)
        ltg_df = runner_ltg_selection.blast_output_to_ltg_tax_id()
        pandas._testing.assert_frame_equal(self.ltg_bak_df, ltg_df)
