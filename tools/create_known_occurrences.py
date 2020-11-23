import os
import pathlib
import shutil
import pandas
import sqlalchemy

from vtam.models.FilterCodonStop import FilterCodonStop
from vtam.models.Marker import Marker
from vtam.models.Run import Run
from vtam.models.Sample import Sample
from vtam.models.SampleInformation import SampleInformation
from vtam.models.Variant import Variant
from vtam.utils.PathManager import PathManager
from sqlalchemy.orm import sessionmaker

## Set variables

test_path = os.path.join(PathManager.get_test_path())
outdir_path = os.path.join(test_path, 'outdir')
shutil.rmtree(outdir_path, ignore_errors=True)
pathlib.Path(outdir_path).mkdir(parents=True, exist_ok=True)

occurrences_keep_tsv = os.path.join(test_path, "test_files_dryad.f40v5_small", "run1_mfzr_zfzr", "known_occurrences_keep.tsv")
sample_type_tsv = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tools", "sample_types.tsv")
sample_information_path = os.path.join(test_path, "test_files_dryad.f40v5_small", "run1_mfzr_zfzr", "sample_information.tsv")
filter_codon_stop_path = os.path.join(test_path, "test_files_dryad.f40v5_small", "run1_mfzr_zfzr", "filter_codon_stop.tsv")
variant_path = os.path.join(test_path, "test_files_dryad.f40v5_small", "run1_mfzr_zfzr", "variant_filter_codon_stop.tsv")
db_sqlite = os.path.join(outdir_path, "db.sqlite")

## Init DB and create session

engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_sqlite), echo=False)

sample_information_df = pandas.read_csv(sample_information_path, sep="\t", header=0)
sample_information_df.to_sql(name=SampleInformation.__tablename__, con=engine.connect(), if_exists='replace')

run_df = pandas.DataFrame({'name': ['run1']}, index=range(1, 2))
run_df.to_sql(name=Run.__tablename__, con=engine.connect(), index_label='id', if_exists='replace')

marker_df = pandas.DataFrame({'name': ['MFZR', 'ZFZR']}, index=range(1, 3))
marker_df.to_sql(name=Marker.__tablename__, con=engine.connect(), index_label='id', if_exists='replace')

sample_df = pandas.DataFrame({'name': ['tpos1_run1', 'tnegtag_run1', '14ben01', '14ben02']}, index=range(1, 5))
sample_df.to_sql(name=Sample.__tablename__, con=engine.connect(), index_label='id', if_exists='replace')

variant_df = pandas.read_csv(variant_path, sep="\t", header=0, index_col='id')
variant_df.to_sql(name=Variant.__tablename__, con=engine.connect(), index_label='id', if_exists='replace')

filter_codon_stop_df = pandas.read_csv(filter_codon_stop_path, sep="\t", header=0)
filter_codon_stop_df.to_sql(name=FilterCodonStop.__tablename__, con=engine.connect(), if_exists='replace')

Session = sessionmaker(bind=engine)
session = Session()

# Update one of the variant as a bad variant
session.query(Customers).filter(Customers.id! = 2).update({Customers.name:"Mr."+Customers.name}, synchronize_session = False)

## Load files
keep_df = pandas.read_csv(occurrences_keep_tsv, sep="\t", header=0)
sample_type_df = pandas.read_csv(sample_type_tsv, sep="\t", header=0)

##

print(session.query(FilterCodonStop).filter_by(sample_id=2).count())

