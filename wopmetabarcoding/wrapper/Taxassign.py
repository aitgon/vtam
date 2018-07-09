from wopmars.framework.database. tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger
from wopmetabarcoding.wrapper.TaxassignUtilities import alignment_vsearch, create_phylogenetic_line_df, sub_fasta_creator,dataframe2ltgdefinition, rank_hierarchy, seq2tax_db_sqlite_to_df, create_tsv_per_variant, get_vsearch_results_per_variant, taxassignation
import pandas,os
from wopmetabarcoding.utils.constants import tempdir
from Bio import SeqIO
from numpy import nan


class Taxassign(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.Taxassign"
    }
    __input_file_taxassign_db = "taxassign_db_fasta"
    __marker_variant_path = "marker_variant_path"
    __assignlvl2id = "assignlvl2id"
    __tax_assign_db_sqlite = "tax_assign_db_sqlite"
    __default_output = 'default_output'
    __udb_database = 'udb_database'

    def specify_input_file(self):
        return [
            Taxassign.__input_file_taxassign_db,
            Taxassign.__marker_variant_path,
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
        marker_variant_path = self.input_file(Taxassign.__marker_variant_path)
        tax_assign_pars_tsv = self.input_file(Taxassign.__assignlvl2id)
        tax_assign_sqlite = self.input_file(Taxassign.__tax_assign_db_sqlite)
        udb_database = self.input_file(Taxassign.__udb_database)
        #
        # 80.0 class order    5
        #
        default_output = self.output_file(Taxassign.__default_output)
        with open(marker_variant_path, 'r') as fin:

            for marker_line in fin:
                marker_line = marker_line.strip().split('\t')
                marker_name = marker_line[0]
                marker_variant_fasta = marker_line[2] # path to fasta with filtered variants
                result_dataframe = pandas.read_csv(marker_line[1], sep="\t")
                result_dataframe["taxa"] = nan
                print(result_dataframe)
                marker_variant_vsearch_tsv = marker_variant_fasta.replace('.fasta', '.tsv')
                # output_tsv = output_tsv.replace(tempdir, '/tmp/tmpe6yiaf0x/')
                nb_variants = 100
                # Loop over groups of records
                sub_fasta_path_list = sub_fasta_creator(marker_variant_fasta, nb_variants, marker_name)
                for sub_fasta_path in sub_fasta_path_list:
                    alignment_vsearch(sub_fasta_path, udb_database, marker_variant_vsearch_tsv)
                    create_tsv_per_variant(marker_variant_vsearch_tsv, "data/dbtaxa.sqlite")
                    for record in SeqIO.parse(sub_fasta_path, 'fasta'):
                        tsv_output = os.path.join(tempdir, (marker_name + "_"  + record.description + '.tsv'))
                        get_vsearch_results_per_variant("data/dbtaxa.sqlite", record.description, tsv_output)
                        taxassignation(tsv_output, tax_assign_sqlite, tax_assign_pars_tsv, result_dataframe, record.description)
                        print(tsv_output)
        print(result_dataframe)








