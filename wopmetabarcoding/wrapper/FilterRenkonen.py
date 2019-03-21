import pandas, itertools


class FilterRankonen:
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.FilterRankonen"
    }

    # Input file
    # Input table
    __input_table_marker = "Marker"
    __input_table_run = "Run"
    __input_table_biosample = "Biosample"
    __input_table_replicate = "Replicate"
    __input_file_sample2fasta = "sample2fasta"
    __input_table_filter_pcr_error = "FilterPCRError"
    __input_table_Variant = "Variant"
    # Output table
    __output_table_filter_chimera = "FilterRankonen"

    def specify_input_file(self):
        return [
            FilterRankonen.__input_file_sample2fasta,

        ]

    def specify_input_table(self):
        return [
            FilterRankonen.__input_table_marker,
            FilterRankonen.__input_table_run,
            FilterRankonen.__input_table_biosample,
            FilterRankonen.__input_table_replicate,
            FilterRankonen.__input_table_filter_pcr_error,
            FilterRankonen.__input_table_Variant,
        ]

    def specify_output_table(self):
        return [
            FilterRankonen.__output_table_filter_chimera,
        ]

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        #
        # Input file path
        input_file_sample2fasta = self.input_file(FilterRankonen.__input_file_sample2fasta)
        #
        # Input table models
        marker_model = self.input_table(FilterRankonen.__input_table_marker)
        run_model = self.input_table(FilterRankonen.__input_table_run)
        pcr_error_model = self.input_table(FilterRankonen.__input_table_filter_pcr_error)
        biosample_model = self.input_table(FilterRankonen.__input_table_biosample)
        replicate_model = self.input_table(FilterRankonen.__input_table_replicate)
        variant_model = self.input_table(FilterRankonen.__input_table_Variant)
        #

        #
        # Output table models
        filter_chimera_model = self.output_table(FilterRankonen.__output_table_filter_chimera)

        ##########################################################
        #
        # 1. Read sample2fasta to get run_id, marker_id, biosample_id, replicate_id for current analysis
        #
        ##########################################################
        sample2fasta_df = pandas.read_csv(input_file_sample2fasta, sep="\t", header=None, \
                                          names=['tag_forward', 'primer_forward', 'tag_reverse', 'primer_reverse',
                                                 'marker_name', 'biosample_name', \
                                                 'replicate_name', 'run_name', 'fastq_fwd', 'fastq_rev', 'fasta'])
        sample_instance_list = []
        for row in sample2fasta_df.itertuples():
            marker_name = row.marker_name
            run_name = row.run_name
            biosample_name = row.biosample_name
            replicate_name = row.replicate_name
            with engine.connect() as conn:
                # get run_id ###########
                stmt_select_run_id = select([run_model.__table__.c.id]).where(run_model.__table__.c.name == run_name)
                run_id = conn.execute(stmt_select_run_id).first()[0]
                # get marker_id ###########
                stmt_select_marker_id = select([marker_model.__table__.c.id]).where(
                    marker_model.__table__.c.name == marker_name)
                marker_id = conn.execute(stmt_select_marker_id).first()[0]
                # get biosample_id ###########
                stmt_select_biosample_id = select([biosample_model.__table__.c.id]).where(
                    biosample_model.__table__.c.name == biosample_name)
                biosample_id = conn.execute(stmt_select_biosample_id).first()[0]
                # get replicate_id ###########
                stmt_select_replicate_id = select([replicate_model.__table__.c.id]).where(
                    replicate_model.__table__.c.name == replicate_name)
                replicate_id = conn.execute(stmt_select_replicate_id).first()[0]
                # add this sample_instance ###########
                sample_instance_list.append({'run_id': run_id, 'marker_id': marker_id, 'biosample_id': biosample_id,
                                             'replicate_id': replicate_id})

        ##########################################################
        #
        # 2. Delete /run/markerbiosample/replicate from this filter table
        #
        ##########################################################
        with engine.connect() as conn:
            conn.execute(filter_chimera_model.__table__.delete(), sample_instance_list)

        ##########################################################
        #
        # 3. Select marker/run/biosample/replicate from variant_read_count_model
        #
        ##########################################################

        pcr_error_model_table = pcr_error_model.__table__
        stmt_variant_filter_lfn = select([pcr_error_model_table.c.marker_id,
                                          pcr_error_model_table.c.run_id,
                                          pcr_error_model_table.c.variant_id,
                                          pcr_error_model_table.c.biosample_id,
                                          pcr_error_model_table.c.replicate_id,
                                          pcr_error_model_table.c.filter_id,
                                          pcr_error_model_table.c.filter_delete,
                                          pcr_error_model_table.c.read_count]) \
            .where(pcr_error_model_table.c.filter_id == 10) \
            .where(pcr_error_model_table.c.filter_delete == 0)
        # Select to DataFrame
        variant_filter_lfn_passed_list = []
        with engine.connect() as conn:
            for row in conn.execute(stmt_variant_filter_lfn).fetchall():
                variant_filter_lfn_passed_list.append(row)
        variant_read_count_df = pandas.DataFrame.from_records(variant_filter_lfn_passed_list,
                                                              columns=['marker_id', 'run_id', 'variant_id',
                                                                       'biosample_id', 'replicate_id', 'filter_id',
                                                                       'filter_delete', 'read_count'])
        # run_id, marker_id, variant_id, biosample_id, replicate_id, read_count, filter_id, filter_delete
        variant_model_table = variant_model.__table__
        stmt_variant = select([variant_model_table.c.id,
                               variant_model_table.c.sequence])

        # Select to DataFrame
        variant_filter_lfn_passed_list = []
        with engine.connect() as conn:
            for row in conn.execute(stmt_variant).fetchall():
                variant_filter_lfn_passed_list.append(row)
        variant_df = pandas.DataFrame.from_records(variant_filter_lfn_passed_list,
                                                   columns=['id', 'sequence'])
        # import pdb;
        # pdb.set_trace()
        ##########################################################
        #
        # 4. Run Filter
        #
        ##########################################################

        # import pdb;
        # pdb.set_trace()
        ##########################################################
        #
        # 5. Insert Filter data
        #
        ##########################################################
        # 1St way to write output on table: error
        # records = df_filter_output.to_dict('records')
        # with engine.connect() as conn:
        #         conn.execute(filter_chimera_model.__table__.insert(), records)
        # second way worked




def  renkonen_distance(variant_read_count_df, run_id, marker_id, biosample_id, left_replicate_id, right_replicate_id):
    #  Compute sum of read_count per 'run_id', 'marker_id', 'biosample_id', 'replicate_id'
    variant_read_proportion_per_replicate_df = variant_read_count_df[
        ['run_id', 'marker_id', 'biosample_id', 'replicate_id', 'read_count']].groupby(
        by=['run_id', 'marker_id', 'biosample_id', 'replicate_id']).sum().reset_index()
    variant_read_proportion_per_replicate_df = variant_read_proportion_per_replicate_df.rename(
        columns={'read_count': 'read_count_sum_per_replicate'})

    # Merge variant read_count with read_count_sum_per_replicate
    variant_read_proportion_per_replicate_df = variant_read_count_df.merge(variant_read_proportion_per_replicate_df,
                                                                           left_on=['run_id', 'marker_id',
                                                                                    'biosample_id', 'replicate_id'],
                                                                           right_on=['run_id', 'marker_id',
                                                                                     'biosample_id', 'replicate_id'])

    variant_read_proportion_per_replicate_df['variant_read_count_propotion_per_replicate'] \
        = variant_read_proportion_per_replicate_df.read_count / variant_read_proportion_per_replicate_df.read_count_sum_per_replicate
    variant_read_proportion_per_replicate_df.drop('read_count', axis=1, inplace=True)
    variant_read_proportion_per_replicate_df.drop('read_count_sum_per_replicate', axis=1, inplace=True)
    #

    # Select the read proportion for the biosample_id, left_replicate
    left_variant_read_proportion_per_replicate_per_biosample_df = variant_read_proportion_per_replicate_df.loc[
        (variant_read_proportion_per_replicate_df.run_id == run_id)
        & (variant_read_proportion_per_replicate_df.marker_id == marker_id)
        & (variant_read_proportion_per_replicate_df.biosample_id == biosample_id)
        & (variant_read_proportion_per_replicate_df.replicate_id == left_replicate_id)
        ]

    # Select the read proportion for the biosample_id, left_replicate
    right_variant_read_proportion_per_replicate_per_biosample_df = variant_read_proportion_per_replicate_df.loc[
        (variant_read_proportion_per_replicate_df.run_id == run_id)
        & (variant_read_proportion_per_replicate_df.marker_id == marker_id)
        & (variant_read_proportion_per_replicate_df.biosample_id == biosample_id)
        & (variant_read_proportion_per_replicate_df.replicate_id == right_replicate_id)
        ]


    #  Merge left and right replicate
    variant_read_proportion_per_replicate_left_right = left_variant_read_proportion_per_replicate_per_biosample_df.merge( \
        right_variant_read_proportion_per_replicate_per_biosample_df,
        on=['run_id', 'marker_id', 'variant_id', 'biosample_id'])
    # rename columns
    variant_read_proportion_per_replicate_left_right = variant_read_proportion_per_replicate_left_right.rename(
        columns={'replicate_id_x': 'replicate_id_left'})
    variant_read_proportion_per_replicate_left_right = variant_read_proportion_per_replicate_left_right.rename(
        columns={'variant_read_count_propotion_per_replicate_x': 'variant_read_count_propotion_per_replicate_left'})
    variant_read_proportion_per_replicate_left_right = variant_read_proportion_per_replicate_left_right.rename(
        columns={'replicate_id_y': 'replicate_id_right'})
    variant_read_proportion_per_replicate_left_right = variant_read_proportion_per_replicate_left_right.rename(
        columns={'variant_read_count_propotion_per_replicate_y': 'variant_read_count_propotion_per_replicate_right'})

    # variant_read_proportion_per_replicate_left_right = variant_read_proportion_per_replicate_left_right[['variant_id','rp_of_variant_in_replicate1', 'rp_of_variant_in_replicate2']]

    variant_read_proportion_per_replicate_left_right['min_read_proportion'] = \
    variant_read_proportion_per_replicate_left_right[
        ['variant_read_count_propotion_per_replicate_left', 'variant_read_count_propotion_per_replicate_right']].apply(
        lambda row: row.min(), axis=1)


    distance_left_right = 1 - sum(variant_read_proportion_per_replicate_left_right['min_read_proportion'])
    return distance_left_right


def renkonen_distance_df(variant_read_count_df, run_id, marker_id, biosample_id, left_replicate_id, right_replicate_id):
    variant_read_proportion_per_replicate_df = variant_read_count_df[
        ['run_id', 'marker_id', 'biosample_id', 'replicate_id', 'read_count']].groupby(
        by=['run_id', 'marker_id', 'biosample_id', 'replicate_id']).sum().reset_index()
    variant_read_proportion_per_replicate_df = variant_read_proportion_per_replicate_df.rename(
        columns={'read_count': 'read_count_sum_per_replicate'})

    # Merge variant read_count with read_count_sum_per_replicate
    variant_read_proportion_per_replicate_df = variant_read_count_df.merge(variant_read_proportion_per_replicate_df,
                                                                           left_on=['run_id', 'marker_id',
                                                                                    'biosample_id', 'replicate_id'],
                                                                           right_on=['run_id', 'marker_id',
                                                                                     'biosample_id', 'replicate_id'])

    variant_read_proportion_per_replicate_df['variant_read_count_propotion_per_replicate'] \
        = variant_read_proportion_per_replicate_df.read_count / variant_read_proportion_per_replicate_df.read_count_sum_per_replicate
    variant_read_proportion_per_replicate_df.drop('read_count', axis=1, inplace=True)
    variant_read_proportion_per_replicate_df.drop('read_count_sum_per_replicate', axis=1, inplace=True)
    #
    run_id = 1
    marker_id = 1
    biosample_id = 1
    left_replicate_id = 1
    right_replicate_id = 2
    # Select the read proportion for the biosample_id, left_replicate
    left_variant_read_proportion_per_replicate_per_biosample_df = variant_read_proportion_per_replicate_df.loc[
        (variant_read_proportion_per_replicate_df.run_id == run_id)
        & (variant_read_proportion_per_replicate_df.marker_id == marker_id)
        & (variant_read_proportion_per_replicate_df.biosample_id == biosample_id)
        & (variant_read_proportion_per_replicate_df.replicate_id == left_replicate_id)
        ]

    # Select the read proportion for the biosample_id, right_replicate
    right_variant_read_proportion_per_replicate_per_biosample_df = variant_read_proportion_per_replicate_df.loc[
        (variant_read_proportion_per_replicate_df.run_id == run_id)
        & (variant_read_proportion_per_replicate_df.marker_id == marker_id)
        & (variant_read_proportion_per_replicate_df.biosample_id == biosample_id)
        & (variant_read_proportion_per_replicate_df.replicate_id == right_replicate_id)
        ]
    #  Merge left and right replicate
    variant_read_proportion_per_replicate_left_right = left_variant_read_proportion_per_replicate_per_biosample_df.merge( \
        right_variant_read_proportion_per_replicate_per_biosample_df,
        on=['run_id', 'marker_id', 'variant_id', 'biosample_id'])
    # rename columns
    variant_read_proportion_per_replicate_left_right = variant_read_proportion_per_replicate_left_right.rename(
        columns={'replicate_id_x': 'replicate_id_left'})
    variant_read_proportion_per_replicate_left_right = variant_read_proportion_per_replicate_left_right.rename(
        columns={'variant_read_count_propotion_per_replicate_x': 'variant_read_count_propotion_per_replicate_left'})
    variant_read_proportion_per_replicate_left_right = variant_read_proportion_per_replicate_left_right.rename(
        columns={'replicate_id_y': 'replicate_id_right'})
    variant_read_proportion_per_replicate_left_right = variant_read_proportion_per_replicate_left_right.rename(
        columns={
            'variant_read_count_propotion_per_replicate_y': 'variant_read_count_propotion_per_replicate_right'})

    # variant_read_proportion_per_replicate_left_right = variant_read_proportion_per_replicate_left_right[['variant_id','rp_of_variant_in_replicate1', 'rp_of_variant_in_replicate2']]

    variant_read_proportion_per_replicate_left_right['min_read_proportion'] = \
        variant_read_proportion_per_replicate_left_right[
            ['variant_read_count_propotion_per_replicate_left',
             'variant_read_count_propotion_per_replicate_right']].apply(
            lambda row: row.min(), axis=1)
    variant_read_proportion_per_replicate_left_right['distance_left_right'] = 1 - sum(
        variant_read_proportion_per_replicate_left_right['min_read_proportion'])
    variant_read_proportion_per_replicate_left_right = \
        variant_read_proportion_per_replicate_left_right[
            ['biosample_id', 'replicate_id_left', 'replicate_id_right', 'distance_left_right']].drop_duplicates()

    return variant_read_proportion_per_replicate_left_right;

