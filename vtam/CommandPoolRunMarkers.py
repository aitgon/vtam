import os
import pandas
import pathlib
import sys

import sqlalchemy
from Bio import SeqIO

from sqlalchemy.ext.automap import automap_base
from vtam.utils.ParamsFile import ParamsFile

from vtam.models.Sample import Sample
from vtam.models.FilterCodonStop import FilterCodonStop
from vtam.utils.AsvTableRunner import AsvTableRunner
from vtam.utils.Logger import Logger
from vtam.utils.NameIdConverter import NameIdConverter
from vtam.utils.PathManager import PathManager
from vtam.utils.RunMarkerFile import RunMarkerFile
from vtam.utils.SequenceClusterer import SequenceClusterer
from vtam.utils.VSearch import VSearch
from vtam.utils.VariantDF import VariantDF
from vtam.utils.VTAMexception import VTAMexception


class CommandPoolRunMarkers(object):
    """Class for the Pool Marker wrapper"""

    def __init__(self, asv_table_df, run_marker_df=None):

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

    def run_vsearch_to_cluster_sequences(self):
        # Define fasta_path tsv_path
        fasta_path = os.path.join(self.tmp_dir, 'variants.fa')
        # Create variant variant_read_count_input_df
        variant_df = self.asv_table_df[[
            'variant_id', 'sequence', 'read_count']].drop_duplicates(inplace=False)
        variant_df.columns = ['id', 'sequence', 'size']
        variant_df.set_index('id', inplace=True)
        # Create fasta_path file from asv_table_df
        variant_df_utils = VariantDF(variant_df)
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
        vsearch_cluster = VSearch(parameters=vsearch_parameters)
        vsearch_cluster.run()
        self.cluster_path = vsearch_output_cluster_path
        return vsearch_output_centroid_fasta, vsearch_output_cluster_path

    def get_vsearch_clusters_to_df(self):
        """
        Analysis vsearch cluster output, which a tsv_path that corresponds to the same tsv_path with ticker 0, 1, 2

        For instance, if self.cluster_path=/tmp/tmpibbwi9oc/test_pool_markers.py/cluster.fa,
        then there are /tmp/tmpibbwi9oc/test_pool_markers.py/cluster.fa0, ...1, ...2, etc with the different clusters

        :return: pandas.DataFrame with columns: variant_id_centroid and variant_id
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

        # Drop some columns
        def are_reads(x):
            return int(sum(x) > 0)

        agg_dic = {}
        for k in ['variant_id', 'run_name', 'marker_name', 'pooled_sequences']:
            agg_dic[k] = lambda x: ','.join(map(str, sorted(list(set(x)))))

        for k in self.sample_names:
            agg_dic[k] = are_reads
        pooled_marker_df = self.cluster_df[[
                                               'centroid_variant_id', 'variant_id', 'run_name',
                                               'marker_name'] + self.sample_names]

        ############################################################################################
        #
        # issue 0002645. Add sequences of cluster members with the exception of the centroid variant
        # Annotate all variants with their sequences
        #
        ############################################################################################

        pooled_marker_df = pooled_marker_df.merge(self.asv_table_df[['variant_id', 'sequence']], left_on='variant_id',
                               right_on='variant_id')
        pooled_marker_df = pooled_marker_df.rename(columns={'sequence': 'pooled_sequences'})

        pooled_marker_df = pooled_marker_df.groupby(
            'centroid_variant_id').agg(agg_dic).reset_index()

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
    def main(cls, db, pooled_marker_tsv, run_marker_tsv, params):

        #######################################################################
        #
        # Parameters
        #
        #######################################################################

        # params_dic = constants.get_params_default_dic()
        params_dic = ParamsFile(params).get_params_dic()

        cluster_identity = params_dic['cluster_identity']

        run_marker_file_obj = RunMarkerFile(tsv_path=run_marker_tsv)

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

        asv_table_runner = AsvTableRunner(variant_read_count_df=variant_read_count_df,
                                          engine=engine, sample_list=sample_list, cluster_identity=cluster_identity)
        asv_table_df = asv_table_runner.create_asvtable_df()
        asv_table_df.rename({'run': 'run_name', 'marker': 'marker_name', 'variant': 'variant_id'}, axis=1, inplace=True)

        ############################################################################################
        #
        # Prefix biosample columns with run name for same biosample name in different runs
        #
        ############################################################################################

        asv_table_2_df = asv_table_df.copy()
        asv_table_2_df_col_list = asv_table_2_df.columns.tolist()

        for biosample_name in asv_table_df.iloc[:, 5:-4].columns.tolist():

            run_name_list = asv_table_df[['run_name', biosample_name]]['run_name'].unique().tolist()

            #  run-biosample not unique
            if len(run_name_list) > 1:

                for run_name in run_name_list:
                    colname = run_name + ',' + biosample_name
                    # asv_table_2_df[run_name + "," + biosample_name] = 0
                    asv_table_i_df = asv_table_df.loc[asv_table_df.run_name == run_name][
                        ['run_name', 'marker_name', 'variant_id', biosample_name]]
                    asv_table_i_df.rename({biosample_name: colname}, axis=1,
                                          inplace=True)
                    asv_table_2_df = asv_table_2_df.merge(asv_table_i_df, on=['run_name', 'marker_name', 'variant_id'], how='left')
                    asv_table_2_df[run_name + ',' + biosample_name].fillna(0, inplace=True)  # fillna
                    asv_table_2_df[run_name + ',' + biosample_name] = asv_table_2_df[colname].astype('int')  # convert to int

                    biosample_name_ix = asv_table_2_df_col_list.index(biosample_name)
                    asv_table_2_df_col_list.insert(biosample_name_ix, colname)

                asv_table_2_df_col_list.remove(biosample_name)

        asv_table_2_df = asv_table_2_df[asv_table_2_df_col_list]
        asv_table_2_df.fillna(0, inplace=True)

        ############################################################################################
        #
        # Reorder columns
        #
        ############################################################################################

        column_list = asv_table_2_df.columns.tolist()
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

        pool_marker_runner = CommandPoolRunMarkers(asv_table_df=asv_table_2_df, run_marker_df=run_marker_df)
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

        #######################################################################
        #
        # To tsv
        #
        #######################################################################

        pooled_marker_df.to_csv(pooled_marker_tsv, sep="\t", index=False)