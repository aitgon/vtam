from wopmars.framework.database. tables.ToolWrapper import ToolWrapper
from wopmetabarcoding.wrapper.TaxassignUtilities import indexed_db_creation

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







