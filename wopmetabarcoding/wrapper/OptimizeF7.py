from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger

# from wopmetabarcoding.wrapper.OptimizeLFNutilities import OptimizeLFNRunner
from sqlalchemy import select
import pandas


class OptimizeF7(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.OptimizeLFN"
    }

    # Input file
    __input_file_positive_variants = "positive_variants"
    # Input table
    __input_table_run = "Run"
    __input_table_marker = "Marker"
    __input_table_biosample = "Biosample"
    __input_table_variant_read_count = "VariantReadCount"
    # Output file
    __output_file_optimize_lfn = "optimize_lfn"


    def specify_input_file(self):
        return[
            OptimizeF7.__input_file_positive_variants,
        ]

    def specify_input_table(self):
        return [
            OptimizeF7.__input_table_marker,
            OptimizeF7.__input_table_run,
            OptimizeF7.__input_table_biosample,
            OptimizeF7.__input_table_variant_read_count,
        ]


    def specify_output_file(self):
        return [
            OptimizeF7.__output_file_optimize_lfn,
        ]

    def specify_params(self):
        return {
        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind

        ##########################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ##########################################################
        #
        # Input file path
        input_file_positive_variants = self.input_file(OptimizeF7.__input_file_positive_variants)
        #
        # Input table models
        run_model = self.input_table(OptimizeF7.__input_table_run)
        marker_model = self.input_table(OptimizeF7.__input_table_marker)
        biosample_model = self.input_table(OptimizeF7.__input_table_biosample)
        variant_read_count_model = self.input_table(OptimizeF7.__input_table_variant_read_count)
        #
        # Output file path
        output_file_optimize_lfn = self.output_file(OptimizeF7.__output_file_optimize_lfn)

        ##########################################################
        #
        # 1. Read input_file_positive_variants to get run_id, marker_id, biosample_id, for current analysis
        #
        ##########################################################
        positive_variant_df = pandas.read_csv(input_file_positive_variants, sep="\t", header=0, \
                                 names=['marker_name', 'run_name', 'biosample_name', 'sample_type',
                                            'variant_id', 'action', 'variant_sequence', 'Note'], index_col=False)

        delete_variant_df = pandas.read_csv(input_file_positive_variants, sep="\t", header=0,\
            names=['marker_name', 'run_name', 'biosample_name', 'variant_id', 'variant_sequence'], index_col=False)

        negative_variant_df = pandas.read_csv(input_file_positive_variants, sep="\t", header=0, \
                             names=['marker_name', 'run_name', 'biosample_name', 'variant_id', 'variant_sequence'],
                             index_col=False)

        #########################################################
            #
            # 1-1 extract only the columns of interests and the variant having delete as action
            #
            ##########################################################

        delete_variant_df = delete_variant_df[delete_variant_df["action"] == "delete"]

        delete_variant_df = delete_variant_df[['marker_name', 'run_name', 'biosample_name', 'variant_id', 'variant_sequence']]

        #add the negative variant to delete_df
        delete_variant_df = delete_variant_df.append(negative_variant_df, ignore_index=True)

        ##########################################################
        #
        # 2. Select run_id, marker_id, variant_id, biosample, replicate where variant, biosample, etc in positive_variants_df
        #  3. Get read_count: N_ijk
        #
        ##########################################################
        variant_read_count_list = []
        with engine.connect() as conn:
            for row in delete_variant_df.itertuples():
                run_name = row.run_name
                marker_name = row.marker_name
                biosample_name = row.biosample_name
                variant_id = row.variant_id
                variant_sequence = row.variant_sequence
                stmt_select = select([
                    run_model.__table__.c.id,
                    marker_model.__table__.c.id,
                    variant_read_count_model.__table__.c.variant_id,
                    biosample_model.__table__.c.id,
                    variant_read_count_model.__table__.c.replicate_id,
                    variant_read_count_model.__table__.c.read_count, ]) \
                    .where(run_model.__table__.c.name == run_name) \
                    .where(variant_read_count_model.__table__.c.run_id == run_model.__table__.c.id) \
                    .where(marker_model.__table__.c.name == marker_name) \
                    .where(variant_read_count_model.__table__.c.marker_id == marker_model.__table__.c.id) \
                    .where(biosample_model.__table__.c.name == biosample_name) \
                    .where(variant_read_count_model.__table__.c.biosample_id == biosample_model.__table__.c.id) \
                    .where(variant_read_count_model.__table__.c.variant_id == variant_id) \
                    .distinct()
                variant_read_count_list = variant_read_count_list + conn.execute(stmt_select).fetchall()
        optimized_lfn_df = pandas.DataFrame.from_records(variant_read_count_list,
                                                         columns=['run_id', 'marker_id', 'variant_id', 'biosample_id',
                                                                  'replicate_id', 'N_ijk'])
        optimized_lfn_df.drop_duplicates(inplace=True)


        ##########################################################
        #
        # 4 - Create data frame containing the variant_id, lfn_read_count_threshold, sample_type (mock1,mock2 ...)
        #
        ##########################################################
        lfn_read_count_threshod_dic = {}
        lfn_read_count_threshod_list = []

        sample_type_list = positive_variant_df.sample_type.tolist()
        for sample_type in sample_type_list:
              df =  positive_variant_df[positive_variant_df["sample_type"] == sample_type]
              df = df[['marker_name', 'run_name', 'biosample_name','sample_type','variant_id', 'variant_sequence']]
              variant_id_list = df.variant_id.tolist()
              for variant_id in variant_id_list:
                 df_by_sample_type_df = optimized_lfn_df[optimized_lfn_df["variant_id"] == variant_id]
                 rdc_list =df_by_sample_type_df["N_ijk"].tolist()
                 rdc_max_min  = [max(rdc_list), min(rdc_list)]
                 if  rdc_list[1] not in rdc_max_min: lfn_read_count_threshold = rdc_list[1]
                 if rdc_list[2] not in rdc_max_min: lfn_read_count_threshold = rdc_list[2]
                 if rdc_list[3] not in rdc_max_min: lfn_read_count_threshold = rdc_list[3]
                 #
                 read_count_number= df_by_sample_type_df["N_ijk"].sum()
                 lfn_read_count_threshod_dic = {}
                 # write in data frame the variant_id with his lfn_read_count_threshold and sample_type
                 lfn_read_count_threshod_dic["variant_id"] = variant_id
                 lfn_read_count_threshod_dic["lfn_read_count_threshold"] = lfn_read_count_threshold
                 lfn_read_count_threshod_dic["sample_type"] = sample_type
                 lfn_read_count_threshod_dic["read_count_number"] = read_count_number
              lfn_read_count_threshod_list.append(lfn_read_count_threshod_dic)
        # data frame containing the ['variant_id','lfn_read_count_threshold','sample_type','read_count_number']
        lfn_read_count_threshod_df = pandas.DataFrame(lfn_read_count_threshod_list)

        #choose the lowest value of the lfn_read_count_threshold from the data frame
        lfn_read_count_threshod_df = lfn_read_count_threshod_df.sort_values('lfn_read_count_threshold', ascending=True)
        lfn_read_count_threshold = lfn_read_count_threshod_df['lfn_read_count_threshold'].head(1).values[0]

        # for the lfn_per_variant_threshold Do the ration  lfn_read_count_threshol / number_read_count and choose alsoo the lowest one

        lfn_read_count_threshod_df["lfn_per_variant_threshold"] = lfn_read_count_threshod_df["lfn_read_count_threshold"]/ lfn_read_count_threshod_df["read_count_number"]
        # sort by lowest value lfn_per_variant_threshold
        lfn_read_count_threshod_df = lfn_read_count_threshod_df.sort_values('lfn_per_variant_threshold', ascending=True)
        #
        lfn_per_variant_threshold = lfn_read_count_threshod_df['lfn_per_variant_threshold'].head(1).values[0]




        ##########################################################
        #
        # 5- Run filter for a series of (lfn_read_count_threshold -lfn_per_variant_threshold) combinations
        #
        ##########################################################


        f7_lfn_delete_absolute_read_count()
        f2_f4_lfn_delete_per_sum_variant()







        ##########################################################
        #
        # 6 -determine a specific lfn_per_variant_threshold specified to the delete variant that steel there despite using the optimal combinaison in 5)-
        #
        ##########################################################



        ##########################################################
        #
        # 7. Write TSV file
        #
        ##########################################################
        # optimized_lfn_df.to_csv(output_file_optimize_lfn, header=True, sep='\t', float_format='%.10f', index=False)
