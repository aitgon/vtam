from wopmars.framework.database. tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger
from wopmetabarcoding.wrapper.TaxassignUtilities import vsearch_command, create_phylogenetic_line_df, sub_fasta_creator,dataframe2ltgdefinition, rank_hierarchy, seq2tax_db_sqlite_to_df, vsearch_output_to_sqlite, get_vsearch_output_for_variant_as_df, taxassignation, indexed_db_creation, otu_tables_creator
import pandas,os
from wopmetabarcoding.utils.constants import tempdir
from Bio import SeqIO
from numpy import nan
import errno

class Taxassign(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.Taxassign"
    }
    __input_file_db_udb = "db_udb"
    __marker_variant_path = "marker_variant_path"
    __assignlvl2id = "assignlvl2id"
    __tax_assign_db_sqlite = "tax_assign_db_sqlite"
    __default_output = 'default_output'
    __otu_table_tsv = 'otu_table_tsv'

    def specify_input_file(self):
        return [
            Taxassign.__input_file_db_udb,
            Taxassign.__marker_variant_path,
            Taxassign.__assignlvl2id,
            Taxassign.__tax_assign_db_sqlite,
        ]

    def specify_output_file(self):
        return [
            Taxassign.__default_output,
            Taxassign.__otu_table_tsv,
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
        db_udb = self.input_file(Taxassign.__input_file_db_udb)
        # Output files
        default_output = self.output_file(Taxassign.__default_output)
        otu_file = self.output_file(Taxassign.__otu_table_tsv)
        # path to the tsv with filtered variants
        marker_variant_path = self.input_file(Taxassign.__marker_variant_path)
        tax_assign_pars_tsv = self.input_file(Taxassign.__assignlvl2id)
        tax_assign_sqlite = self.input_file(Taxassign.__tax_assign_db_sqlite)
        #udb_database = self.option("udb_database")
        #indexed_db_creation(taxassign_db_fasta, udb_database)
        #
        # Parameters
        output_dir_taxassign = self.option("output_dir_taxassign")
        try:
            os.makedirs(output_dir_taxassign)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        #
        # 80.0 class order    5
        #
        default_output = self.output_file(Taxassign.__default_output)
        otu_df = pandas.DataFrame()
        with open(marker_variant_path, 'r') as fin:
            # with open(otu_file, 'w') as fout:
            for marker_line in fin:
                marker_line = marker_line.strip().split('\t')
                marker_name = marker_line[0]
                marker_variant_filter_info_tsv = marker_line[1]
                marker_variant_fasta = marker_line[2] # path to fasta with filtered variants
                marker_variant_filter_info_taxa_df = pandas.read_csv(marker_variant_filter_info_tsv, sep="\t")
                marker_variant_filter_info_taxa_df["taxa"] = nan # add column taxa
                marker_variant_filter_info_taxa_df["marker_name"] = marker_name
                # vsearch output file path
                output_vsearch_marker = os.path.join(output_dir_taxassign, "output_vsearch_{}.tsv".format(marker_name))
                nb_variants = 100 # sequence group for vsearch
                # Loop over groups of records
                sub_fasta_path_list = sub_fasta_creator(marker_variant_fasta, nb_variants, marker_name)
                for sub_fasta_path in sub_fasta_path_list:
                    vsearch_command(sub_fasta_path, db_udb, output_vsearch_marker)
                    # sqlite db path to store vsearch output
                    vsearch_output_variant2taxa_seq2perc_identity_sqlite = os.path.join(tempdir, "vsearch_output_variant2taxa_seq2perc_identity.sqlite")
                    vsearch_output_to_sqlite(output_vsearch_marker, vsearch_output_variant2taxa_seq2perc_identity_sqlite)
                    # retrieve and analyse each variant
                    for variant in SeqIO.parse(sub_fasta_path, 'fasta'):
                        variant_seq = variant.description
                        tsv_output = os.path.join(tempdir, (marker_name + "_"  + variant_seq + '.tsv'))
                        vsearch_output_for_variant_df = get_vsearch_output_for_variant_as_df(vsearch_output_variant2taxa_seq2perc_identity_sqlite, variant_seq)
                        taxassignation(vsearch_output_for_variant_df, tax_assign_sqlite, tax_assign_pars_tsv, marker_variant_filter_info_taxa_df, variant_seq)
                        otu_df = otu_df.append(marker_variant_filter_info_taxa_df, ignore_index=True)
                #marker_variant_filter_info_taxa_df.to_csv(default_output, sep='\t', header=True, index=False)
            otu_tables_creator(otu_df, otu_file)

    # def run(self):
    #     session = self.session()
    #     engine = session._WopMarsSession__session.bind
    #     # Input files
    #     db_udb = self.input_file(Taxassign.__input_file_db_udb)
    #     # Output files
    #     default_output = self.output_file(Taxassign.__default_output)
    #     otu_file = self.output_file(Taxassign.__otu_table_tsv)
    #     # path to the tsv with filtered variants
    #     marker_variant_path = self.input_file(Taxassign.__marker_variant_path)
    #     tax_assign_pars_tsv = self.input_file(Taxassign.__assignlvl2id)
    #     tax_assign_sqlite = self.input_file(Taxassign.__tax_assign_db_sqlite)
    #     #udb_database = self.option("udb_database")
    #     #indexed_db_creation(taxassign_db_fasta, udb_database)
    #     #
    #     # Parameters
    #     output_dir_taxassign = self.option("output_dir_taxassign")
    #     try:
    #         os.makedirs(output_dir_taxassign)
    #     except OSError as exception:
    #         if exception.errno != errno.EEXIST:
    #             raise
    #     #
    #     # 80.0 class order    5
    #     #
    #     default_output = self.output_file(Taxassign.__default_output)
    #
    #     with open(marker_variant_path, 'r') as fin:
    #         with open(otu_file, 'w') as fout:
    #             for marker_line in fin:
    #                 marker_line = marker_line.strip().split('\t')
    #                 marker_name = marker_line[0]
    #                 marker_variant_filter_info_tsv = marker_line[1]
    #                 marker_variant_fasta = marker_line[2] # path to fasta with filtered variants
    #                 marker_variant_filter_info_taxa_df = pandas.read_csv(marker_variant_filter_info_tsv, sep="\t")
    #                 marker_variant_filter_info_taxa_df["taxa"] = nan # add column taxa
    #                 print(marker_variant_filter_info_taxa_df)
    #                 # vsearch output file path
    #                 output_vsearch_marker = os.path.join(output_dir_taxassign, "output_vsearch_{}.tsv".format(marker_name))
    #                 nb_variants = 100 # sequence group for vsearch
    #                 # Loop over groups of records
    #                 sub_fasta_path_list = sub_fasta_creator(marker_variant_fasta, nb_variants, marker_name)
    #                 for sub_fasta_path in sub_fasta_path_list:
    #                     vsearch_command(sub_fasta_path, db_udb, output_vsearch_marker)
    #                     # sqlite db path to store vsearch output
    #                     vsearch_output_variant2taxa_seq2perc_identity_sqlite = os.path.join(tempdir, "vsearch_output_variant2taxa_seq2perc_identity.sqlite")
    #                     vsearch_output_to_sqlite(output_vsearch_marker, vsearch_output_variant2taxa_seq2perc_identity_sqlite)
    #                     # retrieve and analyse each variant
    #                     for variant in SeqIO.parse(sub_fasta_path, 'fasta'):
    #                         variant_seq = variant.description
    #                         tsv_output = os.path.join(tempdir, (marker_name + "_"  + variant_seq + '.tsv'))
    #                         vsearch_output_for_variant_df = get_vsearch_output_for_variant_as_df(vsearch_output_variant2taxa_seq2perc_identity_sqlite, variant_seq)
    #                         taxassignation(vsearch_output_for_variant_df, tax_assign_sqlite, tax_assign_pars_tsv, marker_variant_filter_info_taxa_df, variant_seq)
    #                 marker_variant_filter_info_taxa_df.to_csv(default_output, sep='\t', header=True, index=False)
    #                 otu_tables_creator(marker_variant_filter_info_taxa_df, fout, marker_name)







