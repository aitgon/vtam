import inspect
import os
import pathlib

import pandas
from sqlalchemy import select
from wopmars.models.ToolWrapper import ToolWrapper

from vtam.utils.VariantDFutils import VariantDFutils
from vtam.utils.Logger import Logger
from vtam.utils.PathManager import PathManager
from vtam.utils.TaxAssignRunner import TaxAssignRunner


class TaxAssign(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.TaxAssign"}

    # Input file
    __input_file_fastainfo = "fastainfo"
    __input_file_taxonomy = "taxonomy"
    # __input_file_map_taxids = "map_taxids"
    # Input table
    __input_table_marker = "Marker"
    __input_table_run = "Run"
    __input_table_biosample = "Biosample"
    __input_table_filter_codon_stop = "FilterCodonStop"
    __input_table_variant = "Variant"

    # Output table
    __output_table_tax_assign = "TaxAssign"

    def specify_input_file(self):
        return [
            TaxAssign.__input_file_fastainfo,
            TaxAssign.__input_file_taxonomy,
        ]

    def specify_input_table(self):
        return [
            TaxAssign.__input_table_marker,
            TaxAssign.__input_table_run,
            TaxAssign.__input_table_biosample,
            TaxAssign.__input_table_variant,
            TaxAssign.__input_table_filter_codon_stop,
        ]

    def specify_output_table(self):
        return [
            TaxAssign.__output_table_tax_assign,

        ]

    def specify_params(self):
        return {
            "ltg_rule_threshold": "float",  # percentage
            "include_prop": "float",  # percentage
            "min_number_of_taxa": "int",  # count
            "blast_db": "str",
        }

    def run(self):
        session = self.session
        engine = session._session().get_bind()
        threads = int(os.getenv('VTAM_THREADS'))

        this_temp_dir = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__))
        pathlib.Path(this_temp_dir).mkdir(exist_ok=True)

        #########################################################
        #
        # 1. Wrapper inputs, outputs and parameters
        #
        #########################################################

        Logger.instance().debug(
            "file: {}; line: {}; Wrapper inputs, outputs and parameters.".format(__file__,
                                                                                 inspect.currentframe().f_lineno, ))
        #
        # Input file
        input_file_fastainfo = self.input_file(TaxAssign.__input_file_fastainfo)
        input_file_taxonomy_tsv = self.input_file(TaxAssign.__input_file_taxonomy)
        #
        # Input table models
        marker_model = self.input_table(TaxAssign.__input_table_marker)
        run_model = self.input_table(TaxAssign.__input_table_run)
        biosample_model = self.input_table(TaxAssign.__input_table_biosample)
        filter_codon_stop_model = self.input_table(TaxAssign.__input_table_filter_codon_stop)
        variant_model = self.input_table(TaxAssign.__input_table_variant)
        # Output table models
        tax_assign_model = self.output_table(TaxAssign.__output_table_tax_assign)
        #
        # Options
        ltg_rule_threshold = float(self.option("ltg_rule_threshold"))  # percentage
        include_prop = float(self.option("include_prop"))  # percentage
        min_number_of_taxa = int(self.option("min_number_of_taxa"))  # count
        blast_db_dir = str(self.option("blast_db"))  # count

        ##########################################################
        #
        # 2. Read fastainfo to get run_id, marker_id, biosample_id, replicate for current analysis
        #
        ##########################################################

        fastainfo_df = pandas.read_csv(input_file_fastainfo, sep="\t", header=0, \
                                       names=['tag_forward', 'primer_forward', 'tag_reverse', 'primer_reverse',
                                              'marker_name', 'biosample_name', \
                                              'replicate', 'run_name', 'fastq_fwd', 'fastq_rev', 'fasta_path'])
        sample_instance_list = []
        for row in fastainfo_df.itertuples():
            marker_name = row.marker_name
            run_name = row.run_name
            biosample_name = row.biosample_name
            replicate = row.replicate
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
                # add this sample_instance ###########
                sample_instance_list.append({'run_id': run_id, 'marker_id': marker_id, 'biosample_id': biosample_id,
                                             'replicate': replicate})

        ##########################################################
        #
        # 3a. Select variants from this run/marker/biosample/replicate combination
        # 3b. Delete variants from this run/marker/biosample/replicate combination
        #
        ##########################################################
        #
        # 3a. Select variants from this run/marker/biosample/replicate combination
        #
        variant_id_delete_list = []
        for sample_instance in sample_instance_list:
            stmt_select_variant_id_delete = select([filter_codon_stop_model.__table__.c.variant_id]) \
                .where(filter_codon_stop_model.__table__.c.run_id == sample_instance['run_id']) \
                .where(filter_codon_stop_model.__table__.c.marker_id == sample_instance['marker_id']) \
                .where(filter_codon_stop_model.__table__.c.biosample_id == sample_instance['biosample_id']) \
                .where(filter_codon_stop_model.__table__.c.replicate == sample_instance['replicate'])
            # Select to DataFrame
            with engine.connect() as conn:
                for row in conn.execute(stmt_select_variant_id_delete).fetchall():
                    variant_id = row[0]
                    if not variant_id in variant_id_delete_list:
                        variant_id_delete_list.append(row)
        #
        # 3b. Delete variants from this run/markerbiosample/replicate combination
        #
        variant_instance_list = [{'variant_id': variant_id} for variant_id in variant_id_delete_list]
        with engine.connect() as conn:
            conn.execute(tax_assign_model.__table__.delete(), variant_instance_list)

        ##########################################################
        #
        # Get variants that passed the filter and are not already assigned in TaxAssign
        #
        ##########################################################

        Logger.instance().debug(
            "file: {}; line: {}; Get variants and sequences that passed the filters".format(__file__,
                                                                                            inspect.currentframe().f_lineno,
                                                                                            'TaxAssign'))

        tax_assign_model_table = tax_assign_model.__table__
        stmt_variant_tax_assign = select([tax_assign_model_table.c.variant_id])

        # These are the variants that are already in taxassign and do not need recalculate
        variant_tax_assign_list = []
        with engine.connect() as conn:
            variant_tax_assign_list = [variant_id[0] for variant_id in conn.execute(stmt_variant_tax_assign).fetchall()]

        filter_codon_stop_model_table = filter_codon_stop_model.__table__
        variant_model_table = variant_model.__table__
        stmt_variant = select([filter_codon_stop_model_table.c.variant_id, variant_model_table.c.sequence]) \
            .where(filter_codon_stop_model_table.c.variant_id.notin_(variant_tax_assign_list)) \
            .where(filter_codon_stop_model_table.c.variant_id == variant_model_table.c.id) \
            .where(filter_codon_stop_model_table.c.filter_delete == 0)\
            .distinct()\
            .order_by("variant_id")

        # Select to DataFrame the variants that will be tried to be assigned to taxa
        variant_list = []
        with engine.connect() as conn:
            for row in conn.execute(stmt_variant).fetchall():
                variant_list.append({'id': row.variant_id, 'sequence': row.sequence})
        variant_df = pandas.DataFrame.from_records(variant_list, index='id')

        ##########################################################
        #
        # 2 Create FASTA file with Variants
        #
        ##########################################################

        Logger.instance().debug(
            "file: {}; line: {}; Create Fasta from Variants".format(__file__, inspect.currentframe().f_lineno))
        variant_fasta = os.path.join(this_temp_dir, 'variant.fasta')
        variant_df_utils = VariantDFutils(variant_df)
        variant_df_utils.to_fasta(variant_fasta)

        ##########################################################
        #
        # Run TaxAssignRunner
        #
        ##########################################################

        taxonomy_df = pandas.read_csv(input_file_taxonomy_tsv, sep="\t", header=0,
                                         dtype={'tax_id': 'int', 'parent_tax_id': 'int', 'old_tax_id': 'float'})

        tax_assign_runner = TaxAssignRunner(variant_df=variant_df, taxonomy_df=taxonomy_df,
                                            blast_db_dir=blast_db_dir, ltg_rule_threshold=ltg_rule_threshold,
                                            include_prop=include_prop, min_number_of_taxa=min_number_of_taxa,
                                            num_threads=threads)
        ltg_df = tax_assign_runner.ltg_df

        ##########################################################
        #
        # 6. Insert Filter data
        #
        ##########################################################
        Logger.instance().debug("file: {}; line: {}; Insert variant_id, ltg_tax_id, ltg_rank to DB".format(__file__,
                                                                                                           inspect.currentframe().f_lineno))
        with engine.connect() as conn:
            conn.execute(tax_assign_model.__table__.insert(), ltg_df.to_dict('records'))
