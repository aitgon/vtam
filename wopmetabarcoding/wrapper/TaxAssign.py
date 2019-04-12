import sqlalchemy
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger


from wopmetabarcoding.wrapper.FilterLFNutilities import FilterLFNRunner
from sqlalchemy import select
import pandas


class TaxAssign(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.TaxAssign"
    }

    # Input file
    # Input table
    # Output table


    def specify_input_file(self):
        return[
        ]

    def specify_input_table(self):
        return [
        ]


    def specify_output_table(self):
        return [
        ]

    def specify_params(self):
        return {
        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        #
        # Input file path
        #
        # Input table models
        #
        # Output table models
        #
        # Options
        # TaxAssign parameters
        #
        ##########################################################
        #
        # Create Fasta from Variants that passed the Filters
        # test_f06_1_create_nucl_gb_accession2taxid_sqlite
        # Run qblast: test_f06_2_run_qblast
        # test_f06_3_annotate_blast_output_with_tax_id & test_f03_import_blast_output_into_df
        # test_f03_1_tax_id_to_taxonomy_lineage
        # test_f05_select_ltg_identity
        #
        ##########################################################
