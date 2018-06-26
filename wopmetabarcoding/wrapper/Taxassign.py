from wopmars.framework.database. tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger
from wopmetabarcoding.wrapper.TaxassignUtilities import alignment_vsearch, create_phylogenetic_line_df, sub_fasta_creator,dataframe2ltgdefinition, rank_hierarchy, seq2tax_db_sqlite_to_df, create_tsv_per_variant
import pandas
from Bio import SeqIO


class Taxassign(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.Taxassign"
    }
    __input_file_taxassign_db = "taxassign_db_fasta"
    __filtered_dataframe_path = "filtered_dataframe_path"
    __assignlvl2id = "assignlvl2id"
    __tax_assign_db_sqlite = "tax_assign_db_sqlite"
    __default_output = 'default_output'
    __udb_database = 'udb_database'

    def specify_input_file(self):
        return [
            Taxassign.__input_file_taxassign_db,
            Taxassign.__filtered_dataframe_path,
            Taxassign.__assignlvl2id,
            Taxassign.__tax_assign_db_sqlite,
            Taxassign.__udb_database

        ]

    def specify_output_file(self):
        return [
            Taxassign.__default_output
        ]

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        # Input files
        taxassign_db_fasta = self.input_file(Taxassign.__input_file_taxassign_db)
        # path to the tsv with filtered variants
        filter_output = self.input_file(Taxassign.__filtered_dataframe_path)
        tax_assign_pars_tsv = self.input_file(Taxassign.__assignlvl2id)
        tax_assign_sqlite = self.input_file(Taxassign.__tax_assign_db_sqlite)
        udb_database = self.input_file(Taxassign.__udb_database)
        #
        # 80.0 class order    5
        #
        default_output = self.output_file(Taxassign.__default_output)
        with open(filter_output, 'r') as fin:
            for line in fin:
                line = line.strip().split('\t')
                marker_name = line[0]
                filtered_variants_fasta = line[2] # path to fasta with filtered variants
                output_tsv = filtered_variants_fasta.replace('.fasta', '.tsv')
                # output_tsv = output_tsv.replace(tempdir, '/tmp/tmpe6yiaf0x/')
                nb_variants = 100
                # Loop over groups of records
                file_list = sub_fasta_creator(filtered_variants_fasta, nb_variants, marker_name)
                for file in file_list:
                    # alignment_vsearch(file, udb_database, output_tsv)
                    create_tsv_per_variant(output_tsv, "data/dbtaxa.sqlite")








