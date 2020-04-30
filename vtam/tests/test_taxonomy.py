# -*- coding: utf-8 -*-
import pathlib
import shutil

from vtam import CommandTaxonomy
from vtam.utils.PathManager import PathManager
from unittest import TestCase

import os


class TestTaxonomy(TestCase):

    def setUp(self):

        testdir_path = os.path.join(PathManager.get_test_path())
        self.outdir_path = os.path.join(testdir_path, "outdir")
        pathlib.Path(self.outdir_path).mkdir(exist_ok=True, parents=True)
        self.taxonomy_tsv = os.path.join(self.outdir_path, "taxonomy.tsv")

    def test_precomputed(self):
        CommandTaxonomy(taxonomy_tsv=self.taxonomy_tsv).download_precomputed_taxonomy()
        self.assertTrue(os.path.getsize(self.taxonomy_tsv) >= 115036883)

    # This test is very slow but tests creating denove taxonomy tsv file from NCBI
    # def test_create_denovo_from_ncbi(self):
    #     CommandTaxonomy(taxonomy_tsv=self.taxonomy_tsv).create_denovo_from_ncbi()
    #     self.assertTrue(os.path.getsize(self.taxonomy_tsv) >= 115036883)

    def tearDown(self):
        shutil.rmtree(self.outdir_path, ignore_errors=True)
