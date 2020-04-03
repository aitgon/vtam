import os
import pandas
import pathlib
import sqlalchemy
import sys

from Bio import SeqIO

from sqlalchemy.ext.automap import automap_base
from sqlalchemy import create_engine

from vtam.models.Biosample import Biosample
from vtam.models.FilterChimeraBorderline import FilterChimeraBorderline
from vtam.models.FilterCodonStop import FilterCodonStop
from vtam.models.Marker import Marker
from vtam.models.Run import Run
from vtam.models.Variant import Variant
from vtam.utils.AsvTableRunner import AsvTableRunner
from vtam.utils.Logger import Logger
from vtam.utils.PathManager import PathManager
from vtam.utils.SampleInformationUtils import SampleInformationUtils
from vtam.utils.VSearch import VSearch
from vtam.utils.VariantDFutils import VariantDFutils
from vtam.utils.VTAMexception import VTAMexception


class RunMarkerTSVreader():
    """Prepares different DFs: engine, variant_read_count_input_df, variant_df, run_df, marker_df, biosample_df,
                                          variant_to_chimera_borderline_df based on run_marker_tsv and db"""

    def __init__(self, db, run_marker_tsv):
        self.__db = db
        #
        run_marker_df = pandas.read_csv(run_marker_tsv, sep="\t", header=0)

        engine = create_engine('sqlite:///{}'.format(db), echo=False)
        Base = automap_base()
        Base.prepare(engine, reflect=True)

        run_declarative = Base.classes.Run
        marker_declarative = Base.classes.Marker
        sample_information_declarative = Base.classes.SampleInformation

        conn = engine.connect()

        ################################################################################################################
        #
        # Create SampleInformationUtils object
        #
        ################################################################################################################

        sample_record_list = []
        for row in run_marker_df.itertuples():
            # Get run_id
            run_row = conn.execute(sqlalchemy.select([run_declarative.id]).where(run_declarative.name == row.run)).first()
            if not (run_row is None):
                run_id = run_row[0]
                # Get marker_id
                marker_row = conn.execute(sqlalchemy.select([marker_declarative.id]).where(marker_declarative.name == row.marker)).first()
                if not (marker_row is None):
                    marker_id = marker_row[0]
                    # Get sample_information entries for this run and marker
                    select_query = sqlalchemy.select([sample_information_declarative.run_id,
                                                      sample_information_declarative.marker_id,
                                                      sample_information_declarative.biosample_id,
                                                      sample_information_declarative.replicate,
                                                      ])\
                        .where(sample_information_declarative.run_id == run_id)\
                        .where(sample_information_declarative.marker_id == marker_id).distinct()
                    result_proxy = conn.execute(select_query)
                    for result_proxy_row in result_proxy:
                        sample_record_list.append(dict(zip(result_proxy_row.keys(), result_proxy_row)))
        conn.close()

        sample_information_df = pandas.DataFrame.from_records(data=sample_record_list)
        sample_information_utils = SampleInformationUtils(engine=engine, sample_information_df=sample_information_df)

        ################################################################################################################
        #
        # Compute all variant_read_count_input_df required for ASV table
        #
        ################################################################################################################

        variant_read_count_df = sample_information_utils.get_variant_read_count_df(FilterCodonStop)

        variant_df = sample_information_utils.get_variant_df(FilterCodonStop, Variant)

        run_df = sample_information_utils.get_run_df(Run)

        marker_df = sample_information_utils.get_marker_df(Marker)

        biosample_df = sample_information_utils.get_biosample_df(Biosample)

        variant_to_chimera_borderline_df = sample_information_utils.get_variant_to_chimera_borderline_df(FilterChimeraBorderline)

        asv_table_runner = AsvTableRunner(engine=engine, variant_read_count_df=variant_read_count_df, variant_df=variant_df,
                                          run_df=run_df, marker_df=marker_df, biosample_df=biosample_df,
                                          variant_to_chimera_borderline_df=variant_to_chimera_borderline_df,
                                          taxonomy_tsv=None)
        self.asv_df_final = asv_table_runner.run()


class CommandPoolRunMarkers(object):
    """Class for the Pool Marker wrapper"""

    def __init__(self, asv_table_df, run_marker_df=None):

        try:
            assert asv_table_df.columns.tolist()[:5] == ['variant_id', 'marker', 'run', 'sequence_length',
                                                         'read_count']
            assert asv_table_df.columns.tolist()[-2:] == ['chimera_borderline', 'sequence']
        except:
            Logger.instance().error(VTAMexception("The ASV table structure is wrong. It is expected to start with these columns:"
                                                  "'variant_id', 'marker', 'run', 'sequence_length', 'read_count'"
                                                  " followed by biosample names and ending with"
                                                  "'phylum', 'class', 'order', 'family', 'genus', 'species', 'ltg_tax_id', "
                                                  "'ltg_tax_name', 'identity', 'ltg_rank', 'chimera_borderline', 'sequence'."))
            sys.exit(1)

        self.biosample_names = asv_table_df.columns.tolist()[5:-2]

        if run_marker_df is None:  # Default: pool all marker
            self.asv_table_df = asv_table_df
        else:  # if run_marker_df: pool only markers in this variant_read_count_input_df
            self.asv_table_df = asv_table_df.merge(run_marker_df, on=['run', 'marker'])

        self.tmp_dir = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__))
        pathlib.Path(self.tmp_dir).mkdir(exist_ok=True)

        self.cluster_path = None # returned by run_vsearch_to_cluster_sequences

        self.cluster_df = None # returned by get_vsearch_clusters_to_df

    def run_vsearch_to_cluster_sequences(self):
        # Define fasta_path path
        fasta_path = os.path.join(self.tmp_dir, 'variants.fa')
        # Create variant variant_read_count_input_df
        variant_df = self.asv_table_df[['variant_id', 'sequence', 'read_count']].drop_duplicates(inplace=False)
        variant_df.columns = ['id', 'sequence', 'size']
        variant_df.set_index('id', inplace=True)
        # Create fasta_path file from asv_table_df
        variant_df_utils = VariantDFutils(variant_df)
        variant_df_utils.to_fasta(fasta_path, add_column='size')
        # Define vsearch output path
        vsearch_output_centroid_fasta = os.path.join(self.tmp_dir, 'centroid.fa')
        # Define cluster output path
        vsearch_output_cluster_path = os.path.join(self.tmp_dir, 'cluster.fa')

        #
        # Create object and run vsearch
        vsearch_parameters = {'cluster_size': fasta_path,
                              'clusters':  vsearch_output_cluster_path,
                              'id': 1, 'sizein': None,
                              'centroids': vsearch_output_centroid_fasta,
                              "threads": int(os.getenv('VTAM_THREADS')),
                              }
        vsearch_cluster = VSearch(parameters = vsearch_parameters)
        vsearch_cluster.run()
        self.cluster_path = vsearch_output_cluster_path
        return vsearch_output_centroid_fasta, vsearch_output_cluster_path

    def get_vsearch_clusters_to_df(self):

        """
        Analysis vsearch cluster output, which a path that corresponds to the same path with ticker 0, 1, 2

        For instance, if self.cluster_path=/tmp/tmpibbwi9oc/test_pool_markers.py/cluster.fa,
        then there are /tmp/tmpibbwi9oc/test_pool_markers.py/cluster.fa0, ...1, ...2, etc with the different clusters

        :return: pandas.DataFrame with columns: variant_id_centroid and variant_id
        """
        if self.cluster_path is None:
            self.run_vsearch_to_cluster_sequences()

        clusters_dir = os.path.dirname(self.cluster_path)
        clusters_basename = os.path.basename(self.cluster_path)
        cluster_i = 0
        cluster_i_path = os.path.join(clusters_dir, clusters_basename + str(cluster_i))
        cluster_instance_list = []
        centroid = None
        while pathlib.Path(cluster_i_path).exists():
            seq_j = 0
            for seq_record in SeqIO.parse(cluster_i_path, "fasta"):
                cluster_instance = {}
                variant_id = seq_record.id.split(';')[0]
                if seq_j == 0:
                    centroid = variant_id
                cluster_instance['centroid_variant_id'] = centroid
                cluster_instance['variant_id'] = variant_id
                seq_j += 1
                cluster_instance_list.append(cluster_instance)
            cluster_i += 1
            cluster_i_path = os.path.join(clusters_dir, clusters_basename + str(cluster_i))
        cluster_df = pandas.DataFrame(data=cluster_instance_list)
        cluster_df = cluster_df.astype(int)
        self.cluster_df = cluster_df
        return cluster_df

    def get_pooled_marker_df(self):

        if self.cluster_df is None:
            self.get_vsearch_clusters_to_df()

        # Merge cluster_df with asv_table
        self.cluster_df = self.cluster_df.merge(self.asv_table_df, on='variant_id')
        self.cluster_df.sort_values(by=self.cluster_df.columns.tolist(), inplace=True)
        #
        # Information to keep about each centroid variant
        centroid_df = self.cluster_df.loc[self.cluster_df.centroid_variant_id == self.cluster_df.variant_id,
                                     ['centroid_variant_id']]
        centroid_df.drop_duplicates(inplace=True)
        #
        # Centroid to aggregated variants and biosamples
        self.cluster_df_col = self.cluster_df.columns
        # biosample_columns = list(set(self.cluster_df_col[6:]).difference(set(list(self.cluster_df_col[-12:]))))
        # Drop some columns
        are_reads = lambda x: int(sum(x)>0)
        agg_dic = {}
        for k in ['variant_id', 'run', 'marker']:
            agg_dic[k] = lambda x: ','.join(map(str, sorted(list(set(x)))))
        for k in self.biosample_names:
            agg_dic[k] = are_reads
        pooled_marker_df = self.cluster_df[['centroid_variant_id', 'variant_id', 'run', 'marker'] + self.biosample_names]
        pooled_marker_df = pooled_marker_df.groupby('centroid_variant_id').agg(agg_dic).reset_index()
        pooled_marker_df = pooled_marker_df.merge(centroid_df, on='centroid_variant_id')
        pooled_marker_df = pooled_marker_df.merge(self.asv_table_df[['variant_id', 'sequence']],
                                                        left_on='centroid_variant_id', right_on='variant_id')
        pooled_marker_df.rename({'variant_id_x': 'variant_id'}, axis=1, inplace=True)
        pooled_marker_df.drop(labels = 'variant_id_y', axis=1, inplace=True)
        return pooled_marker_df

    @classmethod
    def main(cls, db, pooled_marker_tsv, run_marker_tsv):
        run_marker_tsv_reader = RunMarkerTSVreader(db=db, run_marker_tsv=run_marker_tsv)
        if not (run_marker_tsv is None):
            run_marker_df = pandas.read_csv(run_marker_tsv, sep="\t", header=0)
        else:
            run_marker_df = None
        pool_marker_runner = CommandPoolRunMarkers(asv_table_df=run_marker_tsv_reader.asv_df_final, run_marker_df=run_marker_df)
        pooled_marker_df = pool_marker_runner.get_pooled_marker_df()
        pooled_marker_df.to_csv(pooled_marker_tsv, sep="\t", index=False)
