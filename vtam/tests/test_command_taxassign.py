import filecmp
import os
import shlex
import shutil
import sys

import sqlalchemy
import subprocess
import unittest

from urllib import request
from vtam.utils.PathManager import PathManager
from vtam.utils.constants import taxonomy_tsv_gz_url, coi_blast_db_gz_url
from wopmars.Base import Base


class TestCommandTaxAssign(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        # vtam needs to be in the tsv_path
        subprocess.run([sys.executable, '-m', 'pip', 'install', '{}/.'.format(PathManager.get_package_path()),
                        '--upgrade'])

        cls.test_path = os.path.join(PathManager.get_test_path())
        cls.outdir_path = os.path.join(cls.test_path, 'outdir')

        cls.args = {}
        cls.args['taxonomy'] = os.path.join(cls.outdir_path, "taxonomy.tsv")
        cls.args['coi_blast_db_dir'] = os.path.join(cls.outdir_path, "coi_blast_db_dir")

        ############################################################################################
        #
        # Run 'vtam taxonomy'
        #
        ############################################################################################

        command = "vtam taxonomy --output {taxonomy} --precomputed".format(**cls.args)
        subprocess.run(shlex.split(command), check=True)

        ############################################################################################
        #
        # Run 'vtam coi_blast_db'
        #
        ############################################################################################

        command = "vtam coi_blast_db --blastdbdir {coi_blast_db_dir} --blastdbname coi_blast_db_20191211 ".format(**cls.args)
        subprocess.run(shlex.split(command), check=True)

    def setUp(self):

        ################################################################################################################
        #
        # VTAM database
        #
        ################################################################################################################

        self.db_path = os.path.join(self.outdir_path, "db.sqlite")
        self.args['db'] = self.db_path
        self.engine = sqlalchemy.create_engine('sqlite:///{}'.format(self.db_path))  # memory db

        Session = sqlalchemy.orm.sessionmaker(bind=self.engine)
        self.session = Session()

        Base.metadata.create_all(self.engine)

        self.asvtable = os.path.join(self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/asvtable_default.tsv")
        self.args['asvtable'] = self.asvtable
        self.asvtable_taxa = os.path.join(self.outdir_path, "asvtable_default_taxa.tsv")
        self.args['asvtable_taxa'] = self.asvtable_taxa
        self.asvtable_taxa_bak = os.path.join(self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/asvtable_default_taxa.tsv")

    def test_01_command_taxonomy(self):

        self.assertTrue(os.path.getsize(self.args['taxonomy']) >= 115036880)
        self.assertTrue(os.path.getsize(self.args['taxonomy']) <= 115036890)

    def test_02_command_coi_blast_db(self):

        self.assertTrue(os.path.getsize(os.path.join(self.args['coi_blast_db_dir'], 'coi_blast_db_20191211.nsi')) >= 536510)
        self.assertTrue(os.path.getsize(os.path.join(self.args['coi_blast_db_dir'], 'coi_blast_db_20191211.nsi')) <= 536520)
        self.assertTrue(set(os.listdir(self.args['coi_blast_db_dir'])) >= {
            'coi_blast_db_20191211.nsd', 'coi_blast_db_20191211.nhr', 'coi_blast_db_20191211.nsq', 'coi_blast_db_20191211.nin', 'coi_blast_db_20191211.nog',
            'coi_blast_db_20191211.nsi'})

    def test_03_taxassign(self):

        # CommandTaxAssign.main(db=self.db_path, mode='unassigned', variants_tsv=self.asvtable,
        #                       output=self.asvtable_taxa, taxonomy_tsv=self.taxonomy_tsv,
        #                       blastdb_dir_path=self.coi_blast_db_dir, blastdbname_str='coi_blast_db')

        ############################################################################################
        #
        # COI Blast DB
        #
        ############################################################################################

        command = "vtam taxassign --variants {asvtable} --output {asvtable_taxa} --db {db} --blastdbdir {coi_blast_db_dir} --blastdbname coi_blast_db_20191211 --taxonomy {taxonomy}".format(**self.args)
        subprocess.run(shlex.split(command), check=True)

        self.assertTrue(filecmp.cmp(self.asvtable_taxa, self.asvtable_taxa_bak, shallow=True))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.outdir_path, ignore_errors=True)
