import os
import pandas
import pathlib
import shlex
import shutil
import sqlalchemy
import subprocess
import sys
import unittest

from vtam.models.FilterCodonStop import FilterCodonStop
from vtam.models.FilterIndel import FilterIndel
from vtam.models.Marker import Marker
from vtam.models.Run import Run
from vtam.models.Sample import Sample
from vtam.models.Variant import Variant
from vtam.utils import pip_install_vtam_for_tests
from vtam.utils.PathManager import PathManager


class TestWrapperFilterCodonStop(unittest.TestCase):

    def setUp(self):

        pip_install_vtam_for_tests()

        self.test_path = PathManager.get_test_path()
        self.package_path = PathManager.get_package_path()
        self.outdir_path = os.path.join(self.test_path, 'outdir')
        shutil.rmtree(self.outdir_path, ignore_errors=True)
        pathlib.Path(self.outdir_path).mkdir(exist_ok=True, parents=True)

        marker_str = "id name\n1 IIICBR"
        run_str = "id name\n1 TAS2"
        sample_str = "id name\n1 S21"
        variant_str = """id sequence
1 ATTGTCAGACACTCCGTACCATTAGGGTGCTGCAGTCGACTAGTCTATTTTAAGCTTACACGTAGCCGGAATTAGTTCATTACTGGGGTCAATTAATATCATAACAACGATCATTAACTAGAGGGCCCCAGGAATGACCTGGGAGAACTTACCGTTATTCGTGTGGGCTGTATTTATTACAGCGTGGTTACTTGTACTGTCTTTACCAGTACTAGCTGGTGCGATTACCATGCTGCTAACAGATAGGAACTAGAATACTAGTTTCTACGACCCGAACGGAGGAGGAGATCCTCTGCTATACCAGCATCTATTC"""
        filter_indel_str = """id run_id marker_id variant_id sample_id replicate read_count filter_delete
1 1 1 1 1 1 50 0"""

        from sqlalchemy import create_engine

        db_path = os.path.join(self.outdir_path, 'db.sqlite')
        self.engine = create_engine('sqlite:///{}'.format(db_path), echo=True)

        from wopmars.Base import Base
        Session = sqlalchemy.orm.sessionmaker(bind=self.engine)
        self.session = Session()

        Base.metadata.create_all(self.engine)

        from io import StringIO

        run_df = pandas.read_csv(StringIO(run_str), sep=" ")
        run_df.to_sql(name=Run.__tablename__, con=self.engine.connect(), if_exists='append', index=False)

        marker_df = pandas.read_csv(StringIO(marker_str), sep=" ")
        marker_df.to_sql(name=Marker.__tablename__, con=self.engine.connect(), if_exists='append', index=False)

        sample_df = pandas.read_csv(StringIO(sample_str), sep=" ")
        sample_df.to_sql(name=Sample.__tablename__, con=self.engine.connect(), if_exists='append', index=False)

        filter_indel_df = pandas.read_csv(StringIO(filter_indel_str), sep=" ")
        filter_indel_df.to_sql(name=FilterIndel.__tablename__, con=self.engine.connect(), if_exists='append', index=False)

        variant_df = pandas.read_csv(StringIO(variant_str), sep=" ")
        variant_df.to_sql(name=Variant.__tablename__, con=self.engine.connect(), if_exists='append', index=False)

        pathlib.Path(os.path.join(self.outdir_path, "params.yml")).touch()

        sortereadinfo_str = """run	marker	sample	replicate	sortedfasta
TAS2	IIICBR	S21	1	TAS2-R1_S1_L001_R1_001_000.fasta"""

        with open(os.path.join(self.outdir_path, "sortedinfo.tsv"), 'w') as fout:
            fout.write(sortereadinfo_str)



    def test01(self):

        cmd = """wopmars tool vtam.wrapper.FilterCodonStop -D sqlite:///db.sqlite -i \"{'file': {'sortedinfo': 'sortedinfo.tsv', 'params': 'params.yml'}, 'table': {'Marker': 'vtam.models.Marker', 'Run': 'vtam.models.Run', 'Sample': 'vtam.models.Sample', 'FilterIndel': 'vtam.models.FilterIndel', 'Variant': 'vtam.models.Variant'}}\" -o \"{'table': {'FilterCodonStop': 'vtam.models.FilterCodonStop'}}\" -P \"{'genetic_code': 5, 'skip_filter_codon_stop': 0}\""""

        if sys.platform.startswith("win"):
            args = cmd
        else:
            args = shlex.split(cmd)
        subprocess.run(args=args, check=True, cwd=self.outdir_path)

        stmt_select = sqlalchemy.select([FilterCodonStop.__table__.c.filter_delete])
        with self.engine.connect() as conn2:
            filter_delete_value = conn2.execute(stmt_select).first()

        self.assertEqual(True, filter_delete_value[0])

    def tearDown(self):
        shutil.rmtree(self.outdir_path, ignore_errors=True)
