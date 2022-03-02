
import pandas

from wopmars.models.ToolWrapper import ToolWrapper
from vtam.utils.FileSampleInformation import FileSampleInformation
from vtam.utils.ModelVariantReadCountLike import ModelVariantReadCountLike


class ReadCountAverageOverReplicates(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.ReadCountAverageOverReplicates"
    }

    # Input file
    # Input table
    __input_table_marker = "Marker"
    __input_table_run = "Run"
    __input_table_sample = "Sample"
    __input_file_sortedinfo = "sortedinfo"
    __input_table_filter_codon_stop = "FilterCodonStop"
    # Output table
    __output_table_filter_consensus = "ReadCountAverageOverReplicates"

    def specify_input_file(self):
        return[
            ReadCountAverageOverReplicates.__input_file_sortedinfo,

        ]

    def specify_input_table(self):
        return [
            ReadCountAverageOverReplicates.__input_table_marker,
            ReadCountAverageOverReplicates.__input_table_run,
            ReadCountAverageOverReplicates.__input_table_sample,
            ReadCountAverageOverReplicates.__input_table_filter_codon_stop,

        ]

    def specify_output_table(self):
        return [
            ReadCountAverageOverReplicates.__output_table_filter_consensus,

        ]

    def specify_params(self):
        return {
        }

    def run(self):
        session = self.session
        engine = session._session().get_bind()
        #
        # Input file output
        fasta_info_tsv = self.input_file(
            ReadCountAverageOverReplicates.__input_file_sortedinfo)
        #
        codon_stop_model = self.input_table(
            ReadCountAverageOverReplicates.__input_table_filter_codon_stop)

        #
        # Output table models
        consensus_model = self.output_table(
            ReadCountAverageOverReplicates.__output_table_filter_consensus)

        # #######################################################################
        # #
        # # 1. Read sortedinfo to get run_id, marker_id, sample_id, replicate for current analysis
        # #
        # #######################################################################
        #
        # # fasta_info_tsv = FastaInformationTSV(engine=engine, fasta_info_tsv=input_file_sortedinfo)
        # sample_info_tsv_obj = FileSampleInformation(tsv_path=input_file_sortedinfo)
        #
        # #######################################################################
        # #
        # # 2. Delete /run_name/markersamples/replicate from this filter table
        # #
        # #######################################################################
        # # with engine.connect() as conn:
        # #     # conn.execute(consensus_model.__table__.delete(), sample_instance_list)
        # #     conn.execute(consensus_model.__table__.delete(), sample_instance_list)
        # #
        # variant_read_count_like_utils = ModelVariantReadCountLike(
        #     variant_read_count_like_model=consensus_model, engine=engine)
        # sample_record_list = sample_info_tsv_obj.to_identifier_df(
        #     engine=engine).to_dict('records')
        # variant_read_count_like_utils.delete_from_db(
        #     sample_record_list=sample_record_list)
        #
        # #######################################################################
        # #
        # # 3. Select marker_name/run_name/sample/replicate from variant_read_count_model
        # #
        # #######################################################################
        #
        # nijk_df = sample_info_tsv_obj.get_nijk_df(
        #     variant_read_count_like_model=codon_stop_model, filter_id=None)
        #
        # # Exit if no variants for analysis
        # try:
        #     assert nijk_df.shape[0] > 0
        # except AssertionError:
        #     sys.stderr.write(
        #         "Error: No variants available for this filter: {}".format(
        #             os.path.basename(__file__)))
        #     sys.exit(1)

        #######################################################################
        #
        # 1. Read sortedinfo to get run_id, marker_id, sample_id, replicate for current analysis
        # 2. Delete marker_name/run_name/sample/replicate from variant_read_count_model
        # 3. Get nijk_df input
        #
        #######################################################################

        sample_info_tsv_obj = FileSampleInformation(tsv_path=fasta_info_tsv)

        sample_info_tsv_obj.delete_from_db(
            engine=engine, variant_read_count_like_model=consensus_model)

        variant_read_count_df = sample_info_tsv_obj.get_nijk_df(
            variant_read_count_like_model=codon_stop_model,
            engine=engine,
            filter_id=None)

        #######################################################################
        #
        # 4. Run Filter
        #
        #######################################################################

        variant_read_count_delete_df = read_count_average_over_replicates(variant_read_count_df)

        #######################################################################
        #
        # Write to DB
        #
        #######################################################################

        record_list = ModelVariantReadCountLike.filter_delete_df_to_dict(
            variant_read_count_delete_df)
        with engine.connect() as conn:

            # Insert new instances
            conn.execute(consensus_model.__table__.insert(), record_list)

        #######################################################################
        #
        # Touch output tables, to update modification date
        #
        #######################################################################

        for output_table_i in self.specify_output_table():
            declarative_meta_i = self.output_table(output_table_i)
            obj = session.query(declarative_meta_i).order_by(
                declarative_meta_i.id.desc()).first()
            session.query(declarative_meta_i).filter_by(
                id=obj.id).update({'id': obj.id})
            session.commit()


def read_count_average_over_replicates(variant_read_count_df):
    """
        Function used to display the read average of the remaining variant
    """
    variants_sequences = variant_read_count_df['variant_id'].copy()
    variants_sequences = list(set(variants_sequences))
    read_average_columns = ['variant', 'read_average']
    read_average_df = pandas.DataFrame(columns=read_average_columns)

    # sum of read_count over variant_id and sample_id
    read_count_sum_over_variant_id_and_sample_id_df = variant_read_count_df.groupby(
        ['run_id', 'marker_id', 'variant_id', 'sample_id']).sum().reset_index()
    read_count_sum_over_variant_id_and_sample_id_df.drop(
        'replicate', axis=1, inplace=True)
    read_count_sum_over_variant_id_and_sample_id_df = read_count_sum_over_variant_id_and_sample_id_df.rename(
        columns={'read_count': 'read_count'})

    #  count of replicate number per variant_id and sample_id
    replicate_count_over_variant_id_and_sample_id_df = variant_read_count_df.groupby(
        ['run_id', 'marker_id', 'variant_id', 'sample_id']).count().reset_index()
    replicate_count_over_variant_id_and_sample_id_df.drop(
        'read_count', axis=1, inplace=True)
    replicate_count_over_variant_id_and_sample_id_df = replicate_count_over_variant_id_and_sample_id_df.rename(
        columns={'replicate': 'replicate_count'})

    # merge
    df_out = read_count_sum_over_variant_id_and_sample_id_df.merge(
        replicate_count_over_variant_id_and_sample_id_df, left_on=(
            'run_id', 'marker_id', 'variant_id', 'sample_id'), right_on=(
            'run_id', 'marker_id', 'variant_id', 'sample_id'))
    df_out['read_count_average'] = df_out.read_count / df_out.replicate_count
    #
    return df_out
