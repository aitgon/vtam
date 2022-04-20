import filecmp
import os
import shutil
import unittest

from vtam import CommandMakeKnownOccurrences
from vtam.utils.PathManager import PathManager

class TestCommandMakeKnownOccurrences(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.test_path = PathManager.get_test_path()
        cls.outdir_path = os.path.join(cls.test_path, 'outdir')

    def setUp(self):

        self.asvTable = os.path.join(self.test_path, "test_files", "asvtable_default_taxa_test.tsv")
        self.sampleTypes = os.path.join(self.test_path, "test_files", "sample_types.tsv")
        self.mockComposition = os.path.join(self.test_path, "test_files", "mock_composition_test.tsv")

        self.habitat_proportion = "0.5"

        self.known_occurrences = os.path.join(self.outdir_path, "known_occurrences.tsv")
        self.missing_occurrences = os.path.join(self.outdir_path, "missing_occurrences.tsv")

    def test_01(self):

        CommandMakeKnownOccurrences.main(asvTable=self.asvTable, sampleTypes=self.sampleTypes, mockComposition=self.mockComposition, known_occurrences=self.known_occurrences, missing_occurrences=self.missing_occurrences, habitat_proportion=self.habitat_proportion)
        self.assertTrue(filecmp.cmpfiles(self.test_path,self.outdir_path,common=['known_occurrences.tsv','missing_occurrences.tsv'],shallow=True))

    def tearDown(self):
        shutil.rmtree(self.outdir_path, ignore_errors=True)
