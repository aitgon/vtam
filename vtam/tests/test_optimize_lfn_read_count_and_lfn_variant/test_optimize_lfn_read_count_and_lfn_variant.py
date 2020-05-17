import pathlib
import shutil
import unittest

import os
import pandas
import sqlalchemy

from vtam.utils.PathManager import PathManager


class TestOptimizeLFNreadCountAndLFNvariant(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # cls.outdir_path = os.path.join(PathManager.get_test_path(), 'outdir')
        # pathlib.Path(cls.outdir_path).mkdir(parents=True, exist_ok=True)
        pass


    def setUp(self):
        # self.variant_read_count_tsv_path = os.path.join(os.path.dirname(__file__), 'variant_read_count.tsv')
        # self.known_occurrences_path = os.path.join(
        #     PathManager.get_package_path(), 'doc/data/known_occurrences.tsv')
        # self.variant_read_count_df = pandas.read_csv(self.variant_read_count_tsv_path, sep="\t", header=0)
        # self.variant_read_count_df['run_id'] = 1
        # self.db = os.path.join(self.outdir_path, "db.sqlite")
        #
        # engine = sqlalchemy.create_engine('sqlite:///{}'.format(self.db), echo=False)
        # # Base = automap_base()
        # # Base.prepare(engine, reflect=True)
        #
        # # import pdb; pdb.set_trace()
        # self.variant_read_count_df.to_sql('VariantReadCount', con=engine.connect(), if_exists='replace')
        #
        # # o = OptimizeLFNreadCountAndLFNvariant(self.known_occurrences_path)
        # # import pdb; pdb.set_trace()
        pass


    def test_01(self):

        pass

    def tearDown(self):
        # shutil.rmtree(self.outdir_path, ignore_errors=True)
        pass
