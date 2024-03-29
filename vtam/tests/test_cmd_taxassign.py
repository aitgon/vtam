import filecmp
import os
import pathlib
import shlex
import shutil
import sys

import sqlalchemy
import subprocess
import unittest

from vtam.utils.Taxonomy import Taxonomy

from vtam.utils import pip_install_vtam_for_tests
from vtam.utils.PathManager import PathManager
from wopmars.Base import Base


class TestCommandTaxAssign(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        # vtam needs to be in the tsv_path
        pip_install_vtam_for_tests()

        cls.test_path = os.path.join(PathManager.get_test_path())
        cls.outdir_path = os.path.join(cls.test_path, 'outdir')

        cls.args = {}
        cls.args['taxonomy'] = os.path.join(cls.outdir_path, "taxonomy.tsv")
        cls.args['coi_blast_db_dir'] = os.path.join(cls.outdir_path, "coi_blast_db_dir")
        pathlib.Path(cls.args['coi_blast_db_dir']).mkdir(exist_ok=True, parents=True)

        ############################################################################################
        #
        # Run 'vtam taxonomy'
        #
        ############################################################################################

        cmd = "vtam taxonomy --output {taxonomy} --precomputed".format(**cls.args)
        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args)

        ############################################################################################
        #
        # Run 'vtam coi_blast_db'
        #
        ############################################################################################

        cmd = "vtam coi_blast_db --blastdbdir {coi_blast_db_dir} --blastdbname coi_blast_db_20200420 ".format(**cls.args)

        # if not (os.path.isfile(os.path.join(cls.args['coi_blast_db_dir'], "coi_blast_db_20200420.nhr"))):
        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args)

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

        self.asvtable = os.path.join(self.test_path, "test_files_dryad.f40v5_small", "run1_mfzr_zfzr", "asvtable_default.tsv")
        self.args['asvtable'] = self.asvtable
        self.asvtable_taxa = os.path.join(self.outdir_path, "asvtable_default_taxa.tsv")
        self.args['asvtable_taxa'] = self.asvtable_taxa
        self.asvtable_taxa_bak = os.path.join(self.test_path, "test_files_dryad.f40v5_small", "run1_mfzr_zfzr", "asvtable_default_taxa.tsv")

    def test_01_command_taxonomy(self):

        self.assertTrue(os.path.getsize(self.args['taxonomy']) >= 115036880)
        self.assertTrue(os.path.getsize(self.args['taxonomy']) <= 115036890)

    def test_02_command_coi_blast_db(self):

        self.assertTrue(os.path.getsize(os.path.join(self.args['coi_blast_db_dir'], 'coi_blast_db_20200420.nsi')) >= 2915100)
        self.assertTrue(os.path.getsize(os.path.join(self.args['coi_blast_db_dir'], 'coi_blast_db_20200420.nsi')) <= 2915200)
        self.assertTrue(set(os.listdir(self.args['coi_blast_db_dir'])) >= {
            'coi_blast_db_20200420.nsd', 'coi_blast_db_20200420.nhr', 'coi_blast_db_20200420.nsq', 'coi_blast_db_20200420.nin', 'coi_blast_db_20200420.nog',
            'coi_blast_db_20200420.nsi'})

    def test_03_taxassign(self):

        cmd = "vtam taxassign --asvtable {asvtable} --output {asvtable_taxa} --db {db} --blastdbdir {coi_blast_db_dir} --blastdbname coi_blast_db_20200420 --taxonomy {taxonomy}".format(**self.args)

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args)

        self.assertTrue(filecmp.cmp(self.asvtable_taxa, self.asvtable_taxa_bak, shallow=True))

    def test_04_taxonomy_wrong_tax_id(self):

        taxonomy = Taxonomy(tsv=self.args['taxonomy'])
        wrong_tax_id=99999999
        taxonomy.get_one_tax_id_lineage(wrong_tax_id)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.outdir_path, ignore_errors=True)
