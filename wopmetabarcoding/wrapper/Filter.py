import errno
from wopmars.utils.Logger import Logger
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger
from wopmetabarcoding.utils.PathFinder import PathFinder

from wopmetabarcoding.utils.constants import tempdir
from wopmetabarcoding.wrapper.FilterLFN import FilterLFNRunner
from sqlalchemy import select
import pandas, os, pickle


class Filter(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.Filter"
    }

    # Input file
    __input_cutoff_file = "file_cutoff"
    __input_genetic_code_file = "genetic_code"
    # Input table
    __input_table_variant = "Variant"
    __input_table_variant_read_count = "VariantReadCount"
    # Output table
    __output_table_variant_filter_lfn = "VariantFilterLFN"


    def specify_input_file(self):
        return[
            Filter.__input_cutoff_file,
            Filter.__input_genetic_code_file
        ]

    def specify_input_table(self):
        return [
            Filter.__input_table_variant,
            Filter.__input_table_variant_read_count,
        ]


    def specify_output_table(self):
        return [
            Filter.__output_table_variant_filter_lfn,
            Filter.__output_table_variant_filter_non_lfn,
        ]

    def specify_params(self):
        return {
            "lfn_per_variant_threshold": "float",
            "lfn_per_replicate_threshold": "float",
            "lfn_per_biosample_per_replicate_threshold": "float",
            "lfn_read_count_threshold": "float",
            "min_replicate_number": "int",
            "pcr_error_var_prop": "float",
            "renkonen_number_of_replicate": "int",
            "renkonen_renkonen_tail": "float"
        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        #
        # Input file path
        cutoff_file_tsv = self.input_file(Filter.__input_cutoff_file)
        #
        # Input table models
        variant_model = self.input_table(Filter.__input_table_variant)
        variant_read_count_model = self.input_table(Filter.__input_table_variant_read_count)
        #
        # Output table models
        variant_filter_lfn_model = self.output_table(Filter.__output_table_variant_filter_lfn)
        variant_filter_non_lfn_model = self.output_table(Filter.__output_table_variant_filter_non_lfn)
        # with engine.connect() as conn:
        #     conn.execute(variant_filter_lfn_model.__table__.delete())
        #
        # Create variant_biosample_replicate_df to run the filters with:
        # variant_id, biosample_id, replicate_id, read_count, variant_sequence
        variant_model_table = variant_model.__table__
        variant_read_count_model_table = variant_read_count_model.__table__
        stmt_marker_id = select([variant_read_count_model_table.c.marker_id]).distinct()
        # Execute filter in a per_marker basis
        with engine.connect() as conn:
            for row in conn.execute(stmt_marker_id).fetchall():
                marker_id = row[0]
                #
                # Create wrapper outdir in a per marker basis
                # tempdir_marker = os.path.join(tempdir, "FilterUtilities", "marker_id_{}".format(marker_id))
                # PathFinder.mkdir_p(tempdir_marker)
                #
                # Create variant_read_count DF for given marker
                stmt_variant_read_count = select([variant_read_count_model_table.c.variant_id,
                                                               variant_read_count_model_table.c.run_id,
                                                               variant_read_count_model_table.c.biosample_id,
                                                               variant_read_count_model_table.c.replicate_id,
                                                               variant_read_count_model_table.c.read_count])\
                    .where(variant_read_count_model_table.c.marker_id == marker_id)
                variant_read_count_list = []
                for row2 in conn.execute(stmt_variant_read_count).fetchall():
                    variant_read_count_list.append(row2)
                variant_read_count_df = pandas.DataFrame.from_records(variant_read_count_list, columns=['variant_id', 'run_id', 'biosample_id', 'replicate_id', 'read_count'])
                #
                # Create variant DF for given marker
                stmt_variant_read_count = select([variant_model_table.c.id, variant_model_table.c.sequence])\
                    .where(variant_read_count_model_table.c.marker_id == marker_id)\
                    .where(variant_read_count_model_table.c.variant_id == variant_model_table.c.id)
                variant_list = []
                for row2 in conn.execute(stmt_variant_read_count).fetchall():
                    variant_list.append(row2)
                variant_df = pandas.DataFrame.from_records(variant_list, columns=['id', 'sequence'])
                #
                lfn_filter_runner = FilterLFNRunner(variant_df, variant_read_count_df, marker_id)
                #
                # Filter parameters
                lfn_per_variant_threshold = self.option("lfn_per_variant_threshold")
                lfn_per_replicate_threshold = self.option("lfn_per_replicate_threshold")
                lfn_per_replicate_series_threshold = self.option("lfn_per_replicate_series_threshold")
                lfn_read_count_threshold = self.option("lfn_read_count_threshold")
                lfn_per_biosample_per_replicate_threshold = self.option("lfn_per_biosample_per_replicate_threshold")
                #
                Logger.instance().info("Launching LFN filter:")
                #
                ############################################
                # Filter 2: f2_f4_lfn_delete_per_sum_variant
                ############################################
                lfn_filter_runner.f2_f4_lfn_delete_per_sum_variant(lfn_per_variant_threshold)
                #
                ############################################
                # Filter  3: f3_f5_lfn_delete_per_sum_variant_replicate
                ############################################
                lfn_filter_runner.f3_f5_lfn_delete_per_sum_variant_replicate(lfn_per_replicate_threshold)

                ############################################
                # Filter 6:  f6_lfn_delete_per_sum_biosample_replicate_delete
                ############################################

                lfn_filter_runner.f6_lfn_delete_per_sum_biosample_replicate(
                    lfn_per_biosample_per_replicate_threshold)

                ############################################
                # Filter  7:f7_lfn_delete_absolute_read_count
                ############################################
                lfn_filter_runner.f7_lfn_delete_absolute_read_count(lfn_read_count_threshold)

                ############################################
                # Filter 8:f8_lfn_delete_do_not_pass_all_filters
                ############################################
                lfn_filter_runner.f8_lfn_delete_do_not_pass_all_filters()

                ############################################
                # Write all LFN Filters
                ############################################
                records = lfn_filter_runner.delete_variant_df.to_dict('records')
                with engine.connect() as conn:
                    conn.execute(variant_filter_lfn_model.__table__.insert(), records)

                ############################################$
                #
                # NON LFN FILTERS
                #
                ############################################

                variants_passed_filters_lfn_df = lfn_filter_runner.delete_variant_df.copy()

                df2=variants_passed_filters_lfn_df.loc[(variants_passed_filters_lfn_df['filter_id'] == 8) & (variants_passed_filters_lfn_df['filter_delete'] == 0)]

                variants_passed_filters_lfn_df=df2
                non_lfn_filter_runner = FilterLFNRunner(variant_df, variants_passed_filters_lfn_df, marker_id)


                # ############################################
                # # Filter 7: Repeatability: f9_delete_min_replicate_number
                # ############################################
                min_replicate_number = 2
                non_lfn_filter_runner.f9_delete_min_replicate_number(min_replicate_number)
                #
                # Write Non LFN Filters
                records = non_lfn_filter_runner.delete_variant_df.to_dict('records')
                with engine.connect() as conn:
                    conn.execute(variant_filter_non_lfn_model.__table__.insert(), records)




                # ###########################################
                # # Filter 9: PCR error
                # ###########################################
                # lfn_filter_runner.f9_pcr_error(self.option("pcr_error_var_prop"), pcr_error_by_sample=True)
                # #
                # ############################################
                # # Filter 10: Chimera
                # ############################################
                # lfn_filter_runner.f10_chimera(chimera_by_sample_replicate=True, engine=engine, variant_model=variant_model)
                # #
                # ############################################
                # # Filter 11: Renkonen
                # ############################################
                # lfn_filter_runner.f11_renkonen(self.option("renkonen_number_of_replicate"), self.option("renkonen_renkonen_tail"))
                # #
                # ############################################
                # # Filter 12: Indel
                # ############################################
                # lfn_filter_runner.f12_indel()
                # #
                # ############################################
                # # Filter: Stop codon
                # ############################################
                # #     df_codon_stop_per_genetic_code = self.codon_stop_dataframe(genetic_code_tsv)
                # #     Logger.instance().info("Launching pseudogene detection with codon stop filter:")
                # #     filter_lfn_runner.codon_stop(df_codon_stop_per_genetic_code, 2, False)
                # #     filter_lfn_runner.consensus()
                # #     variant2sample2replicate2count_df = filter_lfn_runner.filtered_variants()
                # #
                # ################################
                # # Insert into db
                # ################################
                # lfn_filter_runner.passed_variant_df['variant_id'] = lfn_filter_runner.passed_variant_df.index
                # records = list(lfn_filter_runner.passed_variant_df.T.to_dict().values())
                # with engine.connect() as conn:
                #     conn.execute(variant_selected_model.__table__.insert(), records)


