import filecmp
import os
import pandas
import pathlib
import shutil
import sqlalchemy
import unittest

from vtam.models.Sample import Sample
from vtam.models.FilterChimeraBorderline import FilterChimeraBorderline
from vtam.models.Marker import Marker
from vtam.models.Run import Run
from vtam.models.Variant import Variant
from vtam.utils.AsvTableRunner import AsvTableRunner
from vtam.utils.PathManager import PathManager


class TestAsvtableRunner(unittest.TestCase):

    """Will test main commands based on a complete test dataset"""

    def setUp(self):

        self.package_path = os.path.join(PathManager.get_package_path())
        self.test_path = os.path.join(PathManager.get_test_path())
        self.outdir_path = os.path.join(self.test_path, 'outdir')
        # during development of the test, this prevents errors
        shutil.rmtree(self.outdir_path, ignore_errors=True)
        pathlib.Path(self.outdir_path).mkdir(parents=True, exist_ok=True)

        db_path = os.path.join(self.outdir_path, "db.sqlite")
        filter_codon_stop_path = os.path.join(
            self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/filter_codon_stop.tsv")
        variant_path = os.path.join(
            self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/variant_filter_codon_stop.tsv")
        filter_chimera_borderline_path = os.path.join(
            self.test_path, "test_files_dryad.f40v5_small/run1_mfzr_zfzr/filter_chimera_borderline_and_filter_codon_stop.tsv")

        self.engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_path), echo=False)
        run_df = pandas.DataFrame({'name': ['run1']}, index=range(1, 2))
        run_df.to_sql(name=Run.__tablename__, con=self.engine.connect(), index_label='id')

        marker_df = pandas.DataFrame({'name': ['MFZR', 'ZFZR']}, index=range(1, 3))
        marker_df.to_sql(name=Marker.__tablename__, con=self.engine.connect(), index_label='id')

        sample_df = pandas.DataFrame({'name': ['tpos1_run1', 'tnegtag_run1', '14ben01', '14ben02']}, index=range(1, 5))
        sample_df.to_sql(name=Sample.__tablename__, con=self.engine.connect(), index_label='id')

        variant_df = pandas.read_csv(variant_path, sep="\t", header=0, index_col='id')
        variant_df.to_sql(name=Variant.__tablename__, con=self.engine.connect(), index_label='id')

        filter_chimera_borderline_db = pandas.read_csv(filter_chimera_borderline_path, sep="\t", header=0)
        filter_chimera_borderline_db.to_sql(name=FilterChimeraBorderline.__tablename__, con=self.engine.connect())

        self.filter_codon_stop_df = pandas.read_csv(filter_codon_stop_path, sep="\t", header=0)
        self.sample_list = ['tpos1_run1', 'tnegtag_run1', '14ben01', '14ben02']
        self.outdir_path = os.path.join(self.test_path, 'outdir')

    def test_asvtable(self):

        asvtable_default_path = os.path.join(self.outdir_path, 'asvtable_default.tsv')
        asvtable_default_bak_path = os.path.join(os.path.dirname(__file__), "asvtable_default.tsv")
        asvtable_runner = AsvTableRunner(variant_read_count_df=self.filter_codon_stop_df,
                                         engine=self.engine, sample_list=self.sample_list, cluster_identity=0.97)
        asvtable_runner.to_tsv(asvtable_default_path)

        self.assertTrue(filecmp.cmp(asvtable_default_path, asvtable_default_bak_path, shallow=True))

    def tearDown(self):

        shutil.rmtree(self.outdir_path, ignore_errors=True)
