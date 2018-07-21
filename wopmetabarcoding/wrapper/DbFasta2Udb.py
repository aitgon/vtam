from wopmars.framework.database. tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger
from wopmetabarcoding.wrapper.TaxassignUtilities import vsearch_command, create_phylogenetic_line_df, sub_fasta_creator,dataframe2ltgdefinition, rank_hierarchy, seq2tax_db_sqlite_to_df, vsearch_output_to_sqlite, get_vsearch_output_for_variant_as_df, taxassignation, indexed_db_creation, convert_fileinfo_to_otu_df
import pandas,os
from wopmetabarcoding.utils.constants import tempdir
from Bio import SeqIO
from numpy import nan
import errno

class DbFasta2Udb(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.DbFasta2Udb"
    }
    __input_db_fasta = "db_fasta"
    __output_db_udb = "db_udb"

    def specify_input_file(self):
        return [
            DbFasta2Udb.__input_db_fasta,
        ]

    def specify_output_file(self):
        return [
            DbFasta2Udb.__output_db_udb,
        ]

    def specify_params(self):
        return {
            "output_dir_taxassign": "str",
            "udb_database": "str",
        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        # Input files
        db_fasta = self.input_file(DbFasta2Udb.__input_db_fasta)
        # Output files
        db_udb = self.output_file(DbFasta2Udb.__output_db_udb)
        #
        indexed_db_creation(db_fasta, db_udb)







