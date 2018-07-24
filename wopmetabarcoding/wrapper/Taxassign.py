import pickle
import sqlite3

from wopmars.framework.database. tables.ToolWrapper import ToolWrapper
from wopmetabarcoding.wrapper.TaxassignUtilities import vsearch_command, sub_fasta_creator, vsearch_output_to_sqlite, \
    get_vsearch_output_for_variant_as_df, taxassignation, convert_fileinfo_to_otu_df
import pandas,os
from wopmetabarcoding.utils.constants import tempdir
from wopmetabarcoding.utils.logger import logger, LOGGER_LEVEL
from Bio import SeqIO
import errno
from itertools import repeat

from multiprocessing import Pool

import inspect
from wopmetabarcoding.utils.PathFinder import PathFinder

def f_variant2taxid(variant_seq, variant_class, tax_assign_sqlite, tax_assign_pars_tsv):
    marker_name = variant_class[variant_seq][0]
    logger.debug(
        "file: {}; line: {}; marker_name {} variant_seq {}...".format(__file__, inspect.currentframe().f_lineno,
                                                                   marker_name, variant_seq[1:20]))
    vsearch_per_variant_df = variant_class[variant_seq][1]
    tax_id = taxassignation(variant_seq, marker_name, vsearch_per_variant_df, tax_assign_sqlite, tax_assign_pars_tsv)
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
        #
        # Parameters
        output_dir_taxassign = self.option("output_dir_taxassign")
        PathFinder.mkdir_p(output_dir_taxassign)
        marker2filteranalysis2fasta_df = pandas.read_csv(marker_variant_path, sep="\t", names=['marker', 'variantinfo', 'fasta'], index_col=0)
        ####################
        marker_class = {}
        ####################
        for row in marker2filteranalysis2fasta_df.itertuples(index=True, name='Pandas'):
            logger.debug("file: {}; line: {}; row {}".format(__file__, inspect.currentframe().f_lineno, row))
            marker_name = row[0]
            marker_class[marker_name] = {}
            marker_filterinfo_tsv = row[1]
            ####################
            marker_class[marker_name]['marker_filterinfo_tsv'] = marker_filterinfo_tsv
            ####################
            marker_fasta = row[2]
            #
            # split fasta file for vsearch
            fasta_subset_size = 10
            sub_fasta_dir = os.path.join(tempdir, "TaxAssign", marker_name)
            PathFinder.mkdir_p(sub_fasta_dir)
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
            logger.debug("file: {}; line: {}; marker_name {}".format(__file__, inspect.currentframe().f_lineno, marker_name))
            for sub_fasta_path_i,sub_fasta_path in enumerate(marker_class[marker_name]['sub_fasta_path_list']):
                #
                # run vsearch
                vsearch_output_tsv = os.path.join(tempdir, "TaxAssign", marker_name, "vsearch_i_{}.tsv".format(sub_fasta_path_i))
                logger.debug("file: {}; line: {}; sub_fasta_path {}".format(__file__, inspect.currentframe().f_lineno, sub_fasta_path))
                vsearch_command(sub_fasta_path, db_udb, vsearch_output_tsv)
                marker_class[marker_name]['sub_fasta_path_list'][sub_fasta_path]['vsearch_output_tsv'] = vsearch_output_tsv
                #
                # store vsearch output in sqlite db
                vsearch_output_sqlite = os.path.join(tempdir, "TaxAssign", marker_name, "vsearch_i_{}.sqlite".format(sub_fasta_path_i))
                logger.debug("file: {}; line: {}; vsearch_output_sqlite {}".format(__file__, inspect.currentframe().f_lineno,
                                                                                               vsearch_output_sqlite))
                vsearch_output_to_sqlite(vsearch_output_tsv, vsearch_output_sqlite)
                marker_class[marker_name]['sub_fasta_path_list'][sub_fasta_path]['vsearch_output_sqlite'] = vsearch_output_sqlite
        #
        # for each variant retrieve vsearch output as and store variant, marker, and vsearch output df
        variant_class = {}
        for marker_name in marker_class:
            logger.debug("file: {}; line: {}; marker_name {}".format(__file__, inspect.currentframe().f_lineno, marker_name))
            for sub_fasta_path in marker_class[marker_name]['sub_fasta_path_list']:
                for fasta_entry in SeqIO.parse(sub_fasta_path, 'fasta'):
                    variant_seq = fasta_entry.description
                    vsearch_sqlite = marker_class[marker_name]['sub_fasta_path_list'][sub_fasta_path]['vsearch_output_sqlite']
                    vsearch_per_variant_df = get_vsearch_output_for_variant_as_df(vsearch_sqlite, variant_seq)
                    variant_class[variant_seq] = (marker_name, vsearch_per_variant_df)
        #
        # For each variant, carry out parallel taxassignation
        variant_seq_list = sorted(variant_class.keys())
        if LOGGER_LEVEL == 10:
            PathFindermkdir_p(os.path.join(tempdir, "TaxAssign", "129"))
            variant_seq_list_pkl = os.path.join(tempdir, "TaxAssign", "129", "variant_seq_list.pkl")
            with open(variant_seq_list_pkl, 'wb') as f:
                pickle.dump(variant_seq_list, f)
            logger.debug(
                "file: {}; line: {}; Written {}".format(__file__, inspect.currentframe().f_lineno, variant_seq_list_pkl))
            variant_seq_list_txt = os.path.join(tempdir, "TaxAssign", "129", "variant_seq_list.txt")
            with open(variant_seq_list_txt, 'w') as f:
                for item in variant_seq_list:
                    f.write("{}\n".format(item))
            logger.debug(
                "file: {}; line: {}; Written {}".format(__file__, inspect.currentframe().f_lineno, variant_seq_list_txt))
        # Start of Parallel Version: Comment out after debugging
        with Pool() as p:
            variant2marker2taxid_list = p.starmap(f_variant2taxid, zip(variant_seq_list, repeat(variant_class), repeat(tax_assign_sqlite), repeat(tax_assign_pars_tsv)))
        # End of Parallel Version: Comment out after debugging
        # Start of Non Parallel Version: Comment out after debugging
        # variant2marker2taxid_list = []
        # for variant_seq in variant_seq_list:
        #     tax_id = f_variant2taxid(variant_seq, variant_class, tax_assign_sqlite, tax_assign_pars_tsv)
        #     variant2marker2taxid_list.append(tax_id)
        # End of Non Parallel Version: Comment out after debugging
        if LOGGER_LEVEL == 10:
            variant2marker2taxid_list_pkl = os.path.join(tempdir, "TaxAssign", "variant2marker2taxid_list.pkl")
            with open(variant2marker2taxid_list_pkl, 'wb') as f:
                pickle.dump(variant2marker2taxid_list, f)
            logger.debug(
                "file: {}; line: {}; Written {}".format(__file__, inspect.currentframe().f_lineno, variant2marker2taxid_list_pkl))
            variant2marker2taxid_list_txt = os.path.join(tempdir, "TaxAssign", "variant2marker2taxid_list.txt")
            with open(variant2marker2taxid_list_txt, 'w') as f:
                for item in variant2marker2taxid_list:
                    f.write("{}\n".format(item))
            logger.debug(
                "file: {}; line: {}; Written {}".format(__file__, inspect.currentframe().f_lineno, variant2marker2taxid_list_txt))
        variant2marker2taxid_list_df = pandas.DataFrame.from_records(variant2marker2taxid_list, columns=['variant_seq', 'marker', 'tax_id'])
        if LOGGER_LEVEL == 10:
            variant2marker2taxid_list_df_pkl = os.path.join(tempdir, "TaxAssign", "variant2marker2taxid_list_df.pkl")
            variant2marker2taxid_list_df.to_pickle(variant2marker2taxid_list_df_pkl)
            logger.debug(
                "file: {}; line: {}; Written {}".format(__file__, inspect.currentframe().f_lineno, variant2marker2taxid_list_df_pkl))
            variant2marker2taxid_list_df_tsv = os.path.join(tempdir, "TaxAssign", "variant2marker2taxid_list_df.tsv")
            variant2marker2taxid_list_df.to_csv(variant2marker2taxid_list_df_tsv, sep="\t")
            logger.debug(
                "file: {}; line: {}; Written {}".format(__file__, inspect.currentframe().f_lineno, variant2marker2taxid_list_df_tsv))
        #
        # Prepare otu_marker table (Final output table)
        otu_df = pandas.DataFrame()
        for marker_name in marker_class:
            marker_filterinfo_tsv = marker_class[marker_name]['marker_filterinfo_tsv']
            marker_filterinfo_df = pandas.read_csv(marker_filterinfo_tsv, sep="\t")
            otu_marker_df = convert_fileinfo_to_otu_df(marker_filterinfo_df)
            otu_df = pandas.concat([otu_df, otu_marker_df], axis=0, join='outer', sort=False)
        if LOGGER_LEVEL == 10:
            otu_df_pkl = os.path.join(tempdir, "TaxAssign", "otu_df.pkl")
            otu_df.to_pickle(otu_df_pkl)
            logger.debug(
                "file: {}; line: {}; Written {}".format(__file__, inspect.currentframe().f_lineno, otu_df_pkl))
            otu_df_tsv = os.path.join(tempdir, "TaxAssign", "otu_df.tsv")
            otu_df.to_csv(otu_df_tsv, sep="\t")
            logger.debug(
                "file: {}; line: {}; Written {}".format(__file__, inspect.currentframe().f_lineno, otu_df_tsv))
        #
        # Final operations to format otu_table
        # 
        # Convert tax_id to tax_name
        tax_id_list = variant2marker2taxid_list_df.tax_id.tolist()
        con = sqlite3.connect(tax_assign_sqlite)
        cur = con.cursor()
        sql = "SELECT distinct tax_id,tax_name FROM seq2tax2parent WHERE tax_id in ({tax_id_list})".format(
            tax_id_list=",".join(str(tax_id) for tax_id in tax_id_list))
        cur.execute(sql)
        records = cur.fetchall()
        taxid2taxname_dic = {int(k): v for k, v in records}
        cur.close()
        # taxid2taxname_df = pandas.read_sql(con=con, sql=sql).drop_duplicates()
        #
        # Add variant_seq and tax_id to output otu table
        if LOGGER_LEVEL == 10:
            variant_seq_list_pkl = os.path.join(tempdir, "TaxAssign", "variant_seq_list.pkl")
            with open(variant_seq_list_pkl, 'wb') as f:
                pickle.dump(variant_seq_list, f)
            logger.debug(
                "file: {}; line: {}; Written {}".format(__file__, inspect.currentframe().f_lineno, variant_seq_list_pkl))
            variant_seq_list_txt = os.path.join(tempdir, "TaxAssign", "variant_seq_list.txt")
            with open(variant_seq_list_txt, 'w') as f:
                for item in variant_seq_list:
                    f.write("{}\n".format(item))
            logger.debug(
                "file: {}; line: {}; Written {}".format(__file__, inspect.currentframe().f_lineno, variant_seq_list_txt))
        otu_df['marker'] = None
        otu_df['tax_id'] = None
        otu_df['tax_name'] = None
        for row in variant2marker2taxid_list_df.itertuples():
            logger.debug("file: {}; line: {}; row {}".format(__file__, inspect.currentframe().f_lineno, row))
            variant_seq = row.variant_seq
            tax_id = row.tax_id
            otu_df.loc[otu_df.variant_seq == variant_seq, 'marker'] = row.marker
            try:
                tax_name = taxid2taxname_dic[tax_id]
            except KeyError:
                tax_name = None
            otu_df.loc[otu_df.variant_seq == variant_seq, 'tax_id'] = tax_id
            otu_df.loc[otu_df.variant_seq == variant_seq, 'tax_name'] = tax_name
        #
        # move some columns to beginning and write
        cols = otu_df.columns.tolist()
        cols.insert(0, cols.pop(cols.index('variant_seq')))
        cols.insert(0, cols.pop(cols.index('read_average')))
        cols.insert(0, cols.pop(cols.index('tax_name')))
        cols.insert(0, cols.pop(cols.index('tax_id')))
        cols.insert(0, cols.pop(cols.index('marker')))
        otu_df = otu_df.reindex(columns=cols)
        #
        # format, sort and write
        otu_df.fillna(0, inplace=True)
        otu_df = otu_df * 1
        otu_df = otu_df.apply(lambda col: pandas.to_numeric(col, errors='ignore', downcast='unsigned'))
        otu_df.sort_values(by=otu_df.columns.tolist(), inplace=True)
        otu_df.to_csv(otu_file, sep="\t", index=False, header=True)

