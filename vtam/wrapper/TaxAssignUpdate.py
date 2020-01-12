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


class TaxAssignUpdate(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.TaxAssignUpdate"}

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
            TaxAssignUpdate.__input_file_taxonomy,
        ]

    def specify_input_table(self):
        return [
            TaxAssignUpdate.__input_table_variant,
        ]

    def specify_output_table(self):
        return [
            TaxAssignUpdate.__output_table_tax_assign,

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
        input_file_taxonomy_tsv = self.input_file(TaxAssignUpdate.__input_file_taxonomy)
        #
        # Input table models
        variant_model = self.input_table(TaxAssignUpdate.__input_table_variant)
        # Output table models
        tax_assign_model = self.output_table(TaxAssignUpdate.__output_table_tax_assign)
        #
        # Options
        ltg_rule_threshold = float(self.option("ltg_rule_threshold"))  # percentage
        include_prop = float(self.option("include_prop"))  # percentage
        min_number_of_taxa = int(self.option("min_number_of_taxa"))  # count
        blast_db_dir = str(self.option("blast_db"))  # count

        ##########################################################
        #
        # Get and delete variants in tax_assign table
        #
        ##########################################################

        # 3a. Select variants from this run/marker/biosample/replicate combination

        variant_id_and_sequence_list = []
        # for sample_instance in sample_instance_list:
        stmt_select_variant_id_delete = select([tax_assign_model.__table__.c.variant_id,
                                                variant_model.__table__.c.sequence])\
            .where(tax_assign_model.__table__.c.variant_id == variant_model.__table__.c.id)
        # Select to DataFrame
        with engine.connect() as conn:
            for row in conn.execute(stmt_select_variant_id_delete).fetchall():
                variant_id = row[0]
                if not variant_id in variant_id_and_sequence_list:
                    variant_id_and_sequence_list.append(row)
        #
        # 3b. Delete variants from this run/markerbiosample/replicate combination
        #
        variant_instance_list = [{'id': variant_id_and_sequence[0], 'sequence': variant_id_and_sequence[1]}
                                 for variant_id_and_sequence in variant_id_and_sequence_list]
        if len(variant_instance_list) > 0:
            variant_df = pandas.DataFrame.from_records(variant_instance_list, index='id')

            ##########################################################
            #
            # Create FASTA file with Variants
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
            # Delete old and insert new tax assign
            #
            ##########################################################

            Logger.instance().debug("file: {}; line: {}; Delete old and insert new tax assignation"
                                    .format(__file__, inspect.currentframe().f_lineno))
            with engine.connect() as conn:
                conn.execute(tax_assign_model.__table__.delete(), variant_instance_list)
                conn.execute(tax_assign_model.__table__.insert(), ltg_df.to_dict('records'))
