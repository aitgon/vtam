# -*- coding: utf-8 -*-

from unittest import TestCase

import filecmp
import os
import shutil
import sqlalchemy
import unittest
from wopmars.Base import Base

from urllib import request
from vtam.CommandTaxAssign import CommandTaxAssign
from vtam import CommandTaxonomy
from vtam.CommandBlastCOI import CommandBlastCOI
from vtam.utils.PathManager import PathManager
from vtam.models.Variant import Variant
from vtam.models.TaxAssign import TaxAssign
from vtam.utils.constants import taxonomy_tsv_gz_url, coi_blast_db_gz_url


@unittest.skipIf(
    request.urlopen(taxonomy_tsv_gz_url).getcode() != 200,
    "Precomputed taxonomy file not available online!")
@unittest.skipIf(
    request.urlopen(coi_blast_db_gz_url).getcode() != 200,
    "Precomputed COI blast db database not available online!")
class TestCommandTaxAssign(TestCase):

    @classmethod
    def setUpClass(cls):

        cls.outdir_path = os.path.join(PathManager.get_test_path(), 'outdir')

        ################################################################################################################
        #
        # Taxonomy database
        #
        ################################################################################################################

        cls.taxonomy_tsv = os.path.join(cls.outdir_path, "taxonomy.tsv")
        CommandTaxonomy(taxonomy_tsv=cls.taxonomy_tsv).download_precomputed_taxonomy()

        ################################################################################################################
        #
        # COI Blast DB
        #
        ################################################################################################################

        cls.coi_blast_db_dir = os.path.join(cls.outdir_path, "coi_blast_db")
        CommandBlastCOI(coi_blast_db_dir=cls.coi_blast_db_dir).download()

    def setUp(self):

        ################################################################################################################
        #
        # VTAM database
        #
        ################################################################################################################

        self.db_path = os.path.join(self.outdir_path, "db.sqlite")
        self.engine = sqlalchemy.create_engine('sqlite:///{}'.format(self.db_path))  # memory db

        Session = sqlalchemy.orm.sessionmaker(bind=self.engine)
        self.session = Session()

        Base.metadata.create_all(self.engine)

        self.asvtable = os.path.join(os.path.dirname(__file__), "asvtable.tsv")
        self.asvtable_taxa = os.path.join(self.outdir_path, "asvtable_taxa.tsv")
        self.asvtable_taxa_bak = os.path.join(os.path.dirname(__file__), "asvtable_taxa.tsv")

    def test_01_command_taxonomy(self):

        self.assertTrue(os.path.getsize(self.taxonomy_tsv) >= 115036880)
        self.assertTrue(os.path.getsize(self.taxonomy_tsv) <= 115036890)

    def test_02_command_coi_blast_db(self):

        self.assertTrue(os.path.getsize(os.path.join(self.coi_blast_db_dir, 'coi_blast_db.nsi')) >= 536510)
        self.assertTrue(os.path.getsize(os.path.join(self.coi_blast_db_dir, 'coi_blast_db.nsi')) <= 536520)
        self.assertTrue(set(os.listdir(self.coi_blast_db_dir)) >= {
            'coi_blast_db.nsd', 'coi_blast_db.nhr', 'coi_blast_db.nsq', 'coi_blast_db.nin', 'coi_blast_db.nog',
            'coi_blast_db.nsi'})

    def test_03_taxassign(self):

        CommandTaxAssign.main(db=self.db_path, mode='unassigned', variants_tsv=self.asvtable,
                              output=self.asvtable_taxa, taxonomy_tsv=self.taxonomy_tsv,
                              blasdb_dir_path=self.coi_blast_db_dir, blastdbname_str='coi_blast_db')

        self.assertTrue(filecmp.cmp(self.asvtable_taxa, self.asvtable_taxa_bak, shallow=True))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.outdir_path, ignore_errors=True)
