import os
import pandas
import pathlib
import sys

import sqlalchemy
from Bio import SeqIO

from sqlalchemy.ext.automap import automap_base
from vtam.utils.FileParams import FileParams

from vtam.models.Sample import Sample
from vtam.models.SampleInformation import SampleInformation
from vtam.models.Run import Run
from vtam.models.FilterCodonStop import FilterCodonStop
from vtam.utils.RunnerAsvTable import RunnerAsvTable
from vtam.utils.Logger import Logger
from vtam.utils.NameIdConverter import NameIdConverter
from vtam.utils.PathManager import PathManager
from vtam.utils.FileRunMarker import FileRunMarker
from vtam.utils.SequenceClusterer import SequenceClusterer
from vtam.utils.RunnerVSearch import RunnerVSearch
from vtam.utils.DataframeVariant import DataframeVariant
from vtam.utils.VTAMexception import VTAMexception


class CommandPoolRunMarkers(object):

    """Runner class for the 'pool' command"""

    def __init__(self, asv_table_df, readcounts, run_marker_df=None):
        """
        Constructor of the CommandPoolRunMarkers class

        Parameters
        ----------
        asv_table_df : pandas dataframe
            ASV table.
        readcount : bool
            Default false.
            If false, boolean 0/1 is given for presence or absence of variant in pooled table.
            If true, read integer is given with sum or reads in the pooled runs or markers.
        run_marker_df: pandas dataframe
            Output ASV table with pooled variants
        """

        header = {'run_name', 'marker_name', 'variant_id', 'sequence_length', 'read_count'}
        if not set(asv_table_df.columns) >= header:  # contains at least the 'header_lower' columns
            Logger.instance().error(
                VTAMexception(
                    "The ASV table structure is wrong. It is expected to contain these columns: "
                    "run_name, marker_name, variant_id, sequence_length, read_count"))
            sys.exit(1)

        self.sample_names = asv_table_df.columns.tolist()[5:-2]

        if run_marker_df is None:  # Default: pool all marker_name
            self.asv_table_df = asv_table_df
        else:  # if run_marker_df: pool only markers in this variant_read_count_input_df
            self.asv_table_df = asv_table_df.merge(
                run_marker_df, on=['run_name', 'marker_name'])

        self.tmp_dir = os.path.join(
            PathManager.instance().get_tempdir(),
            os.path.basename(__file__))
        pathlib.Path(self.tmp_dir).mkdir(exist_ok=True)

        self.cluster_path = None  # returned by run_vsearch_to_cluster_sequences

        self.cluster_df = None  # returned by get_vsearch_clusters_to_df
        self.readcounts = readcounts  # returned by get_vsearch_clusters_to_df

    def run_vsearch_to_cluster_sequences(self):
        # Define fasta_path tsv_path
        fasta_path = os.path.join(self.tmp_dir, 'variants.fa')
        # Create variant variant_read_count_input_df
        variant_df = self.asv_table_df[[
            'variant_id', 'sequence', 'read_count']].drop_duplicates(inplace=False)
        variant_df.columns = ['id', 'sequence', 'size']
        variant_df.set_index('id', inplace=True)
        # Create fasta_path file from asv_table_df
        variant_df_utils = DataframeVariant(variant_df)
        variant_df_utils.to_fasta(fasta_path, add_column='size')
        # Define vsearch output tsv_path
        vsearch_output_centroid_fasta = os.path.join(
            self.tmp_dir, 'centroid.fa')
        # Define cluster output tsv_path
        vsearch_output_cluster_path = os.path.join(self.tmp_dir, 'cluster.fa')
        #
        # Create object and run_name vsearch
        vsearch_parameters = {'cluster_size': fasta_path,
                              'clusters': vsearch_output_cluster_path,
                              'id': 1, 'sizein': None,
                              'centroids': vsearch_output_centroid_fasta,
                              "threads": int(os.getenv('VTAM_THREADS')),
                              }
        vsearch_cluster = RunnerVSearch(parameters=vsearch_parameters)
        vsearch_cluster.run()
        self.cluster_path = vsearch_output_cluster_path
        return vsearch_output_centroid_fasta, vsearch_output_cluster_path

    def get_vsearch_clusters_to_df(self):
        """Analysis vsearch cluster output, which a tsv_path that corresponds to the same tsv_path with ticker 0, 1, 2

        For instance, if self.cluster_path=/tmp/tmpibbwi9oc/test_pool_markers.py/cluster.fa,
        then there are /tmp/tmpibbwi9oc/test_pool_markers.py/cluster.fa0, ...1, ...2, etc with the different clusters

        :return:

        Returns
        -------
        pandas.DataFrame with columns: variant_id_centroid and variant_id
        """

        if self.cluster_path is None:
            self.run_vsearch_to_cluster_sequences()

        clusters_dir = os.path.dirname(self.cluster_path)
        clusters_basename = os.path.basename(self.cluster_path)
        cluster_i = 0
        cluster_i_path = os.path.join(
            clusters_dir, clusters_basename + str(cluster_i))
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
            cluster_i_path = os.path.join(
                clusters_dir, clusters_basename + str(cluster_i))
        cluster_df = pandas.DataFrame(data=cluster_instance_list)
        cluster_df = cluster_df.astype(int)
        self.cluster_df = cluster_df
        return cluster_df

    def get_pooled_marker_df(self):

        """Prepares output df"""

        if self.cluster_df is None:
            self.get_vsearch_clusters_to_df()

        # Merge cluster_df with asv_table
        self.cluster_df = self.cluster_df.merge(
            self.asv_table_df, on='variant_id')
        self.cluster_df.sort_values(
            by=self.cluster_df.columns.tolist(), inplace=True)
        #
        #  Information to keep about each centroid variant
        centroid_df = self.cluster_df.loc[self.cluster_df.centroid_variant_id ==
                                          self.cluster_df.variant_id, ['centroid_variant_id']]
        centroid_df.drop_duplicates(inplace=True)

        pooled_marker_df = self.cluster_df[['centroid_variant_id', 'variant_id', 'run_name', 'marker_name'] + self.sample_names]

        ############################################################################################
        #
        # issue 0002645. Add sequences of cluster members with the exception of the centroid variant
        # Annotate all variants with their sequences
        #
        ############################################################################################

        pooled_marker_df = pooled_marker_df.merge(self.asv_table_df[['run_name', 'marker_name', 'variant_id', 'sequence']],
                                                  on=['run_name', 'marker_name', 'variant_id'])
        pooled_marker_df.fillna(0, inplace=True)
        pooled_marker_df = pooled_marker_df.rename(columns={'sequence': 'pooled_sequences'})

        ############################################################################################
        #
        # Dictionnary with aggregation function
        # 'variant_id', 'run_name', 'marker_name', 'pooled_sequences': concatenante unique values
        # biosamples: aggregation with presence:absence or read_sum according to 'readcounts' options
        #
        ############################################################################################

        # Drop some columns
        def are_reads(x):
            if self.readcounts:  # outputs sum
                return int(sum(x))
            else:  # outputs absence/presence
                return int(sum(x) > 0)

        agg_dic = {}
        for k in ['variant_id', 'run_name', 'marker_name', 'pooled_sequences']:
            agg_dic[k] = lambda x: ','.join(map(str, sorted(list(set(x)))))

        for k in self.sample_names:
            agg_dic[k] = are_reads

        # Aggregate by centroid_variant_id

        pooled_marker_df = pooled_marker_df.groupby('centroid_variant_id').agg(agg_dic).reset_index()

        pooled_marker_df = pooled_marker_df.merge(centroid_df, on='centroid_variant_id')
        pooled_marker_df = pooled_marker_df.merge(self.asv_table_df[['variant_id', 'sequence']],
                                                  left_on='centroid_variant_id',
                                                  right_on='variant_id')
        pooled_marker_df.drop_duplicates(inplace=True)
        pooled_marker_df.rename(
            {'variant_id_x': 'variant_id'}, axis=1, inplace=True)
        pooled_marker_df.drop(labels='variant_id_y', axis=1, inplace=True)
        pooled_marker_df.rename(
            {'run_name': 'run', 'marker_name': 'marker', 'centroid_variant_id': 'variant',
             'variant_id': 'pooled_variants'}, axis=1, inplace=True)

        # Move pooled_sequences to end
        sequence_pool_list = pooled_marker_df.pooled_sequences.tolist()
        pooled_marker_df.drop('pooled_sequences', axis=1, inplace=True)
        pooled_marker_df['pooled_sequences'] = sequence_pool_list
        # Move (centroid) sequence to end
        sequence_centroid_list = pooled_marker_df.sequence.tolist()
        pooled_marker_df.drop('sequence', axis=1, inplace=True)
        pooled_marker_df['sequence'] = sequence_centroid_list
        return pooled_marker_df

    @classmethod
    def main(cls, db, pooled_marker_tsv, run_marker_tsv, params, readcounts):

        """

        Parameters
        ----------
        db : str
            path to DB in sqlite format
        pooled_marker_tsv : str
            path to output pooled_marker.tsv file
        run_marker_tsv : str
            path to input run_marker.tsv file
        params
        readcounts: bool
            Output absence/presence (False) or sum of read counts (True)

        Returns
        -------

        """

        #######################################################################
        #
        # Parameters
        #
        #######################################################################

        # params_dic = constants.get_params_default_dic()
        params_dic = FileParams(params).get_params_dic()

        cluster_identity = params_dic['cluster_identity']

        run_marker_file_obj = FileRunMarker(tsv_path=run_marker_tsv)

        # run_marker_tsv_reader = RunMarkerTSVreader(db=db, run_marker_tsv_path=run_marker_tsv)
        if not (run_marker_tsv is None):
            run_marker_df = run_marker_file_obj.read_tsv_into_df()
        else:
            run_marker_df = None

        engine = sqlalchemy.create_engine(
            'sqlite:///{}'.format(db), echo=False)
        Base = automap_base()
        Base.prepare(engine, reflect=True)

        sample_list = run_marker_file_obj.get_sample_ids(engine)
        sample_list = NameIdConverter(id_name_or_sequence_list=sample_list, engine=engine).to_names(Sample)

        ############################################################################################
        #
        # Compute all variant_read_count_input_df required for ASV table
        #
        ############################################################################################

        variant_read_count_df = run_marker_file_obj.get_variant_read_count_df(
            engine=engine, variant_read_count_like_model=FilterCodonStop)

        asv_table_runner = RunnerAsvTable(variant_read_count_df=variant_read_count_df,
                                          engine=engine, sample_list=sample_list, cluster_identity=cluster_identity)
        asv_table_df = asv_table_runner.create_asvtable_df()
        asv_table_df.rename({'run': 'run_name', 'marker': 'marker_name', 'variant': 'variant_id'}, axis=1, inplace=True)

        ############################################################################################
        #
        # Prefix biosample columns with run name for same biosample name in different runs
        #
        ############################################################################################

        asv_table_2_df = asv_table_df.copy()

        for run_name_i, run_name in enumerate(asv_table_df.run_name.unique()):
            asv_table_runi_df = (asv_table_df.loc[asv_table_df.run_name == run_name]).copy()

            for biosample in asv_table_runi_df.iloc[:, 5:-4].columns.tolist():
                asv_table_runi_df.rename({biosample: run_name + '-' + biosample}, axis=1, inplace=True)

            if run_name_i == 0:
                asv_table_2_df = asv_table_runi_df
            else:

                asv_table_2_df = pandas.concat([asv_table_2_df, pandas.DataFrame(columns=asv_table_runi_df.columns)])
                asv_table_2_df = asv_table_2_df.fillna(0)
                asv_table_2_df = pandas.concat([asv_table_2_df, asv_table_runi_df], axis=0, join='outer')

            del (asv_table_runi_df)

        ############################################################################################
        #
        # Reorder columns
        #
        ############################################################################################

        column_list = asv_table_2_df.columns.tolist()
        column_list.remove("run_name")
        column_list.insert(0, "run_name")

        column_list.remove("clusterid")
        column_list.remove("clustersize")
        column_list.remove("chimera_borderline")
        column_list.remove("sequence")
        column_list = column_list + ['clusterid', 'clustersize', 'chimera_borderline', 'sequence']

        column_list.remove("sequence_length")
        column_list.remove("read_count")
        column_list.insert(3, "sequence_length")
        column_list.insert(4, "read_count")

        asv_table_2_df = asv_table_2_df[column_list]

        ############################################################################################
        #
        # Pool markers
        #
        ############################################################################################

        pool_marker_runner = CommandPoolRunMarkers(asv_table_df=asv_table_2_df, run_marker_df=run_marker_df, readcounts=readcounts)
        pooled_marker_df = pool_marker_runner.get_pooled_marker_df()

        #######################################################################
        #
        # Cluster sequences
        #
        #######################################################################

        # reset asvtable-based clusterid and clustersize
        pooled_marker_df.drop(['clusterid', 'clustersize'], axis=1, inplace=True)
        pooled_marker_df.rename({'variant': 'variant_id'}, axis=1, inplace=True)  # prepare
        pooled_marker_df['read_count'] = pooled_marker_df.iloc[:, 4:-2].sum(axis=1)  # prepare

        seq_clusterer_obj = SequenceClusterer(pooled_marker_df, cluster_identity=cluster_identity)
        cluster_count_df = seq_clusterer_obj.compute_clusters()

        pooled_marker_df = pooled_marker_df.merge(cluster_count_df, on='variant_id')
        pooled_marker_df.drop(['read_count'], axis=1, inplace=True)
        ############################################################################################
        #
        # Reorder columns
        #
        ############################################################################################

        column_list = pooled_marker_df.columns.tolist()
        column_list.remove("pooled_sequences")
        column_list.remove("sequence")
        column_list = column_list + ['pooled_sequences', 'sequence']
        pooled_marker_df = pooled_marker_df[column_list]

        # change dtypes
        for col in pooled_marker_df.columns[4:-4]:
            pooled_marker_df[col] = pooled_marker_df[col].astype(int)

        # verify here if the run-sample exists in the sampleinformation database
        # and drop if not
        run_biosample_cols = pooled_marker_df.columns[4:-4]
        # run_biosample_item = run_biosample_cols[0]
        from sqlalchemy.orm import sessionmaker
        Session = sessionmaker(bind=engine)
        session = Session()
        for run_biosample_item in run_biosample_cols:
            thisrun, thisbiosample = run_biosample_item.split('-')
            rowcount = session.query(SampleInformation, Sample, Run).filter(
                SampleInformation.sample_id == Sample.id).filter(
                SampleInformation.run_id == Run.id).filter(
                Run.name == thisrun).filter(Sample.name == thisbiosample).count()
            if rowcount <= 0:
                pooled_marker_df.drop([run_biosample_item], axis=1,
                                      inplace=True)

        #######################################################################
        #
        # To tsv
        #
        #######################################################################

        pooled_marker_df.to_csv(pooled_marker_tsv, sep="\t", index=False)
