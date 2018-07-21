import time
from wopmars.framework.database. tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger
from wopmetabarcoding.wrapper.TaxassignUtilities import vsearch_command, create_phylogenetic_line_df, sub_fasta_creator,dataframe2ltgdefinition, rank_hierarchy, seq2tax_db_sqlite_to_df, vsearch_output_to_sqlite, get_vsearch_output_for_variant_as_df, taxassignation, indexed_db_creation, convert_fileinfo_to_otu_df
import pandas,os
from wopmetabarcoding.utils.constants import tempdir
from Bio import SeqIO
from numpy import nan
import errno
from itertools import repeat

from multiprocessing import Pool

def mkdir_p(path):
    """ Does not fail if directory already exists"""
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def f_variant2taxid(variant_seq, variant_class, marker_class, tax_assign_sqlite, tax_assign_pars_tsv):
    marker_name = variant_class[variant_seq][0]
    vsearch_per_variant_df = variant_class[variant_seq][1]
    # marker_filterinfo_tsv = marker_class[marker_name]['marker_filterinfo_tsv']
    tax_id = taxassignation(variant_seq, vsearch_per_variant_df, tax_assign_sqlite, tax_assign_pars_tsv)
    return (variant_seq, marker_name, tax_id)


class Taxassign(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.Taxassign"
    }
    __input_file_db_udb = "db_udb"
    __marker_variant_path = "marker_variant_path"
    __assignlvl2id = "assignlvl2id"
    __tax_assign_db_sqlite = "tax_assign_db_sqlite"
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
        mkdir_p(output_dir_taxassign)
        #
        otu_table_long_df = pandas.DataFrame()
        marker2filteranalysis2fasta_df = pandas.read_csv(marker_variant_path, sep="\t", names=['marker', 'variantinfo', 'fasta'], index_col=0)
        ####################
        marker_class = {}
        ####################
        for row in marker2filteranalysis2fasta_df.itertuples(index=True, name='Pandas'):
            marker_name = row[0]
            marker_class[marker_name] = {}
            marker_filterinfo_tsv = row[1]
            ####################
            marker_class[marker_name]['marker_filterinfo_tsv'] = marker_filterinfo_tsv
            ####################
            marker_fasta = row[2]
            #
            # split fasta file for vsearch
            fasta_subset_size = 100
            sub_fasta_dir = os.path.join(tempdir, "TaxAssign", marker_name)
            mkdir_p(sub_fasta_dir)
            # sub_fasta_path_template = "marker_%s_i_{}.fasta"%marker_name
            sub_fasta_path_list = sub_fasta_creator(marker_fasta, fasta_subset_size, sub_fasta_dir)
            ####################
            marker_class[marker_name]['sub_fasta_path_list'] = {}
            for sub_fasta_path in sub_fasta_path_list:
                marker_class[marker_name]['sub_fasta_path_list'][sub_fasta_path] = {}
            ####################
        #
        # Run vsearch and and store it in sqlite DB
        for marker_name in marker_class:
            marker_class[marker_name]
            for sub_fasta_path_i,sub_fasta_path in enumerate(marker_class[marker_name]['sub_fasta_path_list']):
                #
                # run vsearch
                vsearch_output_tsv = os.path.join(tempdir, "TaxAssign", marker_name, "vsearch_i_{}.tsv".format(sub_fasta_path_i))
                vsearch_command(sub_fasta_path, db_udb, vsearch_output_tsv)
                marker_class[marker_name]['sub_fasta_path_list'][sub_fasta_path]['vsearch_output_tsv'] = vsearch_output_tsv
                #
                # store vsearch output in sqlite db
                vsearch_output_sqlite = os.path.join(tempdir, "TaxAssign", marker_name, "vsearch_i_{}.sqlite".format(sub_fasta_path_i))
                vsearch_output_to_sqlite(vsearch_output_tsv, vsearch_output_sqlite)
                marker_class[marker_name]['sub_fasta_path_list'][sub_fasta_path]['vsearch_output_sqlite'] = vsearch_output_sqlite
        #
        # for each variant retrieve vsearch output as and store variant, marker, and vsearch output df
        variant_class = {}
        for marker_name in marker_class:
            for sub_fasta_path in marker_class[marker_name]['sub_fasta_path_list']:
                for fasta_entry in SeqIO.parse(sub_fasta_path, 'fasta'):
                    variant_seq = fasta_entry.description
                    vsearch_sqlite = marker_class[marker_name]['sub_fasta_path_list'][sub_fasta_path]['vsearch_output_sqlite']
                    vsearch_per_variant_df = get_vsearch_output_for_variant_as_df(vsearch_sqlite, variant_seq)
                    variant_class[variant_seq] = (marker_name, vsearch_per_variant_df)
        #
        # For each variant, carry out parallel taxassignation
        variant_seq_list = variant_class.keys()
        with Pool() as p:
            # variant_class, marker_class, tax_assign_sqlite, tax_assign_pars_tsv
            variant2marker2taxid_list = p.starmap(f_variant2taxid, zip(variant_seq_list, repeat(variant_class), repeat(marker_class), repeat(tax_assign_sqlite), repeat(tax_assign_pars_tsv)))
        variant2marker2taxid_list_df = pandas.DataFrame.from_records(variant2marker2taxid_list, columns=['variant_seq', 'marker', 'tax_id'])
        #
        # merge tax_ids to other informations and write
        for variant_seq in variant_seq_list:
            marker_name = variant_class[variant_seq][0]
            # vsearch_per_variant_df = variant_class[variant_seq][1]
            marker_filterinfo_tsv = marker_class[marker_name]['marker_filterinfo_tsv']
        #     # tax_id = taxassignation(variant_seq, vsearch_per_variant_df, tax_assign_sqlite, tax_assign_pars_tsv)
        #     #######################""
        #     tax_id = f_variant2taxid(variant_seq, variant_class, marker_class, tax_assign_sqlite, tax_assign_pars_tsv)
        #     print(tax_id)
        #
            marker_filterinfo_df = pandas.read_csv(marker_filterinfo_tsv, sep="\t")
            variant_otu_df = convert_fileinfo_to_otu_df(marker_filterinfo_df)
            if otu_table_long_df.shape == (0,0):
                otu_table_long_df = variant_otu_df
            else:
                # vertical merge
                otu_table_long_df = pandas.concat([otu_table_long_df, variant_otu_df], axis=0, join='outer', sort=False)
        otu_table_long_df = pandas.merge(variant2marker2taxid_list_df, otu_table_long_df, on=['variant_seq'], how='outer')
        otu_table_long_df = otu_table_long_df.reindex(columns=sorted(otu_table_long_df.columns))
        #
        # move some columns to beginning and write
        cols = otu_table_long_df.columns.tolist()
        cols.insert(0, cols.pop(cols.index('variant_seq')))
        cols.insert(0, cols.pop(cols.index('read_average')))
        cols.insert(0, cols.pop(cols.index('tax_id')))
        cols.insert(0, cols.pop(cols.index('marker')))
        otu_table_long_df = otu_table_long_df.reindex(columns=cols)
        #
        # sort and write
        otu_table_long_df.sort_values(by=otu_table_long_df.columns.tolist(), inplace=True)
        otu_table_long_df.to_csv(otu_file, sep="\t", index=False, header=True)
        #
        #
        #
