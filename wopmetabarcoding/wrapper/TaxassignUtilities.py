import inspect

import pandas, os, sqlite3, itertools, csv
from Bio import SeqIO
from numpy import nan
from wopmars.utils.Logger import Logger

from wopmetabarcoding.utils.constants import tempdir
from wopmetabarcoding.utils.logger import logger, LOGGER_LEVEL
from wopmetabarcoding.utils.PathFinder import PathFinder

rank_hierarchy =['no rank', 'phylum', 'superclass', 'class', 'subclass', 'infraclass', 'superorder', 'order', 'suborder', 'infraorder', 'family', 'subfamily', 'genus', 'subgenus', 'species', 'subspecies']


# def create_info_df(fasta_name):
#     columns_name = ['target', 'tax_id', 'name', 'rank', 'parent_id']
#     df_info = pandas.DataFrame(columns=columns_name)
#     for record in SeqIO.parse(fasta_name, 'fasta'):
#         sequence_id = record.description.strip().split('=')
#         seq_name = sequence_id[0].replace(' name', '')
#         name = sequence_id[1].replace(' tax_id', '')
#         tax_id = sequence_id[2].replace(' rank', '')
#         rank = sequence_id[3].replace(' parent_taxid', '')
#         parent_id = sequence_id[4]
#         sequence_info = [seq_name, name, tax_id, rank, parent_id]
#         df_info.loc[len(df_info)] = sequence_info
#     return df_info

# test1 = open(output_tsv, 'w')
                # test1.close()
                # vsearch_usearch_global_args = {'db': taxassign_db_fasta,
                #                                'usearch_global': filtered_variants_fasta,
                #                                'id': str(0.80),
                #                                'maxrejects': 0,
                #                                'maxaccepts': 0,
                #                                'userout': output_tsv,
                #                                'userfields': "--userfields query+target+id --id",
                #                                }
                # vsearch_1 = VSearch1(**vsearch_usearch_global_args)
                # vsearch_1.run()

def indexed_db_creation(taxassign_db_fasta, udb_database):
    """
    Function creating an udb database from the fasta database file for Vsearch to avoid the long indexation step
    at every vsearch
    :param taxassign_db_fasta: Fasta file corresponding at our database for vsearch
    :param udb_database: Fasta file converted to ubd database
    :return: void
    """
    os.system(
        "vsearch --makeudb_usearch " + taxassign_db_fasta + " --output " + udb_database
    )


def vsearch_command(filtered_variants_fasta, taxassign_db_fasta, output_tsv):
    """
    Function charged of the vsearch alignment
    :param filtered_variants_fasta: Fasta file with the sequences of variant which passed the filters.
    :param taxassign_db_fasta:
    :param output_tsv:
    :return:
    """
    os.system(
        "vsearch --usearch_global " + filtered_variants_fasta + " --db " + taxassign_db_fasta +
        " --maxaccept 0 --maxreject 0  --userout " + output_tsv + " --userfields query+target+id --id "
        + str(0.8)
    )


def most_common(lst):
    return max(set(lst), key=lst.count)


def create_list_of_lists(n, tax_seq_id_list):
    tax_seq_id_list = iter(tax_seq_id_list)
    return list(iter(lambda: list(itertools.islice(tax_seq_id_list, n)), []))


def seq2tax_db_sqlite_to_df(seq2tax_db_sqlite, tax_seq_id_list):
    """
    :param seq2tax_db_sqlite:
    :param tax_seq_id_list:
    :return:
    """
    conn = sqlite3.connect(seq2tax_db_sqlite)
    cur = conn.cursor()
    divided_tax_seq_id_list = create_list_of_lists(100, tax_seq_id_list)
    seq2tax_dic_list = []
    for tax_sublist in divided_tax_seq_id_list:
        sql = "select tax_seq_id, tax_name, tax_id, rank_name, tax_parent_id from seq2tax2parent where tax_seq_id in ({seq})".format(
            seq=','.join(['?'] * len(tax_sublist)))
        cur.execute(sql, tax_sublist)
        for row in cur:
            seq2tax_dic_list.append(row)
    seq2tax_df = pandas.DataFrame.from_records(seq2tax_dic_list, columns=["tax_seq_id", "tax_name", "tax_id", "rank_name", "tax_parent_id"])
    cur.close()
    conn.close()
    return seq2tax_df

def create_phylogenetic_line_df(tax_seq_id_list, tax_assign_sqlite):
    """
    Given a list taxon sequence ids (tax_seq_id_list), create a df with a taxon lineage per df

    :param tax_seq_id_list: Integer with taxon sequence id, eg: [7524480]
    :param tax_assign_sqlite: SQLITE DB with tax_seq_id, tax_id and parent_id
    :return: Df with one taxon lineage df per row, eg:
tax_seq_id	no rank	phylum	class	subclass	infraclass	order	suborder	infraorder	family	subfamily	genus	species
5345503	131567	6656	50557	7496	33340	7088	41191	41196	7128	82617	119289	325869
    """
    # Do not log tax_seq_id_list because this is a huge list
    # logger.debug("file: {}; line: {}; tax_seq_id_list {}".format(__file__, inspect.currentframe().f_lineno, tax_seq_id_list))
    conn = sqlite3.connect(tax_assign_sqlite)
    # conn2 = sqlite3.connect(metabarcoding_sqlite)
    lineage_list = []
    tax_lineage_header = ['tax_seq_id'] + rank_hierarchy
    for tax_seq_id in tax_seq_id_list:
        tax_lineage = dict(zip(tax_lineage_header, [None]*len(tax_lineage_header)))
        tax_lineage['tax_seq_id'] = tax_seq_id
        cur = conn.cursor()
        sql = "SELECT tax_id, tax_name, rank_name, tax_parent_id FROM seq2tax2parent WHERE tax_seq_id = ?"
        filter_string = tax_seq_id
        cur.execute(sql, (filter_string,))
        row = cur.fetchone()
        tax_id = row[0]
        rank_name = row[2]
        tax_parent_id = row[3]
        # while not row is None or row[0] != 1:
        while row[0] != 1:
            tax_lineage[rank_name] = tax_id
            filter_string = tax_parent_id
            sql = "SELECT tax_id, tax_name, parent_id, rank FROM tax2parent WHERE tax_id = ?"
            cur.execute(sql, (filter_string,))
            row = cur.fetchone()
            if row is None:
                break
            tax_id = row[0]
            rank_name = row[3]
            tax_parent_id = row[2]
        lineage_list.append(tax_lineage)
        cur.close()
    conn.close()
    tax_lineage_df = pandas.DataFrame(lineage_list)
    tax_lineage_df = tax_lineage_df[tax_lineage_header]
    tax_lineage_df.fillna(value=nan, inplace=True)
    tax_lineage_df = tax_lineage_df.dropna(axis=1, how='all')
    #
    if LOGGER_LEVEL == 10:
        PathFinder.mkdir_p(os.path.join(tempdir, "TaxassignUtilities"))
        tax_lineage_df_pkl = os.path.join(tempdir,"TaxassignUtilities", "tax_lineage_df.pkl")
        tax_lineage_df.to_pickle(tax_lineage_df_pkl)
        logger.debug(
            "file: {}; line: {}; Written {}".format(__file__, inspect.currentframe().f_lineno, tax_lineage_df_pkl))
        tax_lineage_df_tsv = os.path.join(tempdir, "TaxassignUtilities", "tax_lineage_df.tsv")
        tax_lineage_df.to_csv(tax_lineage_df_tsv, sep="\t")
        logger.debug(
            "file: {}; line: {}; Written {}".format(__file__, inspect.currentframe().f_lineno, tax_lineage_df_tsv))
    return tax_lineage_df


def dataframe2ltgdefinition(tax_lineage_df):
    tax_count_perc = pandas.DataFrame({'tax_id': tax_lineage_df.apply(lambda x: x.value_counts().index[0], axis=0)})
    tax_count_perc['count'] = tax_lineage_df.apply(lambda x: x.value_counts().iloc[0], axis=0)
    tax_count_perc['perc'] = tax_count_perc['count'] / tax_lineage_df.shape[0] * 100
    tax_count_perc.drop(['tax_seq_id', 'no rank'], axis=0, inplace=True)
    tax_count_perc = tax_count_perc.loc[tax_count_perc.perc >= 90.0]
    return tax_count_perc


def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.__next__()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


def sub_fasta_creator(fasta_file_path, fasta_subset_size, sub_fasta_dir):
    """
    Split a FASTA file into FASTA pieces of given sequence number

    :param fasta_file_path: Path to FASTA file with marker variants
    :param fasta_subset_size: Maximal sequence number in split FASTA files
    :param marker_name: Marker name to prefix split FASTA files
    :return: List with paths to split FASTA files
    """
    record_iter = SeqIO.parse(open(fasta_file_path), "fasta")
    # print(str(record_iter))
    sub_fasta_path_list = []
    for i, batch in enumerate(batch_iterator(record_iter, fasta_subset_size)):
        sub_fasta_path = os.path.join(sub_fasta_dir, "i_{}.fasta".format(i))
        # filename = filename + "_" + marker_name + ".fasta"
        # sub_fasta_path = os.path.join(tempdir, sub_fasta_path)
        with open(sub_fasta_path, "w") as handle:
            count = SeqIO.write(batch, handle, "fasta")
        sub_fasta_path_list.append(sub_fasta_path)
    # print(sub_fasta_path_list)
    return sub_fasta_path_list

def vsearch_output_to_sqlite(filename, db_to_create):
    """
    Function creating a sqlite table to store data from Vsearch and allow us to filter variant by variant for taxonomic associaton
    :param filename: Vsearch output file
    :param db_to_create:
    :return:
    """
    conn = sqlite3.connect(db_to_create)
    conn.execute("DROP TABLE IF EXISTS alignedtsv")
    cur = conn.cursor()
    conn.execute(
        "CREATE TABLE alignedtsv (id INTEGER PRIMARY KEY AUTOINCREMENT, query_variant VARCHAR, target VARCHAR, identity_thresold FLOAT)"
    )
    with open(filename, 'r') as fin:
        lines = fin.readlines()
        liness = [i.strip().split('\t') for i in lines]
        list_of_lines = create_list_of_lists(10000, liness)
        for lists in list_of_lines:
            cur.executemany("INSERT INTO alignedtsv (query_variant, target, identity_thresold) VALUES (?,?,?);", lists)
    cur.close()
    conn.commit()
    conn.close()


def get_vsearch_output_for_variant_as_df(db_sqlite, variant_seq):
    con = sqlite3.connect(db_sqlite)
    # cur = conn.cursor()
    # data = cur.execute("SELECT query_variant, target, identity_thresold FROM alignedtsv WHERE query_variant == ?", (record_name,))
    # with open(output_tsv, 'w', newline="") as fout:
    #     writer = csv.writer(fout, delimiter='\t')
    #     writer.writerows(data)
    # cur.close()
    sql = "SELECT target, identity_thresold FROM alignedtsv WHERE query_variant == '{}'".format(variant_seq)
    vsearch_output_for_variant_df = pandas.read_sql(sql=sql, con=con)
    con.close()
    return vsearch_output_for_variant_df


def taxassignation(variant_seq, marker_name, vsearch_output_for_variant_df, tax_assign_sqlite, tax_assign_pars_tsv):
    ltg_tax_id = 0 # default
    # vsearch2seq2tax_df = pandas.read_csv(output_tsv, sep="\t", header=None, index_col=1)
    vsearch_output_for_variant_df.columns = ["tax_seq_id", "alignment_identity"]
    # vsearch2seq2tax_df.index.name = 'tax_seq_id'
    #
    tax_seq_id_list = vsearch_output_for_variant_df.tax_seq_id.tolist()
    #
    seq2tax_df = seq2tax_db_sqlite_to_df(tax_assign_sqlite, tax_seq_id_list)
    #
    # tax_assign_pars df
    names = ["identity_threshold", "min_tax_level", "max_tax_resolution", "min_tax_n"]
    tax_assign_pars_df = pandas.read_csv(tax_assign_pars_tsv, sep="\t", header=None, names=names)
    #
    # Merge of the vsearch alignment, the sequence and taxa information
    vsearch_output_for_variant_df[["tax_seq_id"]] = vsearch_output_for_variant_df[["tax_seq_id"]].astype('int64')
    vsearch_output_for_variant_df = pandas.merge(vsearch_output_for_variant_df, seq2tax_df, left_on="tax_seq_id",
                                      right_on="tax_seq_id")
    vsearch_output_for_variant_df = vsearch_output_for_variant_df.assign(
        rank_id=vsearch_output_for_variant_df.rank_name.apply(lambda x: rank_hierarchy.index(x)))
    #
    # Loop over each identity threshold
    for tax_assign_pars_df_row_i, tax_assign_pars_df_row in tax_assign_pars_df.iterrows():
        # identity_threshold, min_tax_level, max_tax_resolution, min_tax_n = 100, "species", "subspecies", 1
        # identity_threshold, min_tax_level, max_tax_resolution, min_tax_n = 97, "genus", "species", 1
        # identity_threshold, min_tax_level, max_tax_resolution, min_tax_n = 95, "family", "species", 3
        # identity_threshold, min_tax_level, max_tax_resolution, min_tax_n = 90, "order", "family", 3
        # identity_threshold, min_tax_level, max_tax_resolution, min_tax_n = 85, "order", "order", 3
        # identity_threshold, min_tax_level, max_tax_resolution, min_tax_n = 80, "class", "order", 5
        identity_threshold, min_tax_level, max_tax_resolution, min_tax_n = tax_assign_pars_df_row.tolist()
        # Logger.instance().info("Selecting sequences with " + str(identity_threshold) + "% identity.")
        min_tax_level_id = rank_hierarchy.index(min_tax_level)
        max_tax_resolution_id = rank_hierarchy.index(max_tax_resolution)
        #
        # test identity_threshold
        # import pdb; pdb.set_trace()
        vsearch2seq2tax_df_selected = vsearch_output_for_variant_df.loc[
            vsearch_output_for_variant_df.alignment_identity >= identity_threshold]
        logger.debug(
            "file: {}; line: {}; marker_name {} variant_seq {}... identity_threshold {} vsearch2seq2tax_df_selected shape {}".format(__file__, inspect.currentframe().f_lineno,
            marker_name, variant_seq[1:20], identity_threshold, vsearch2seq2tax_df_selected.shape))
        if vsearch2seq2tax_df_selected.empty:  #  no lines selected at this alignment identity threshold
            continue  #  next identity threshold
        logger.debug(
            "file: {}; line: {}; marker_name {} variant_seq {}... identity_threshold {} passed".format(__file__, inspect.currentframe().f_lineno,
                                                                       marker_name, variant_seq[1:20], identity_threshold))
        #  continue only if selected lines
        #
        # test min_tax_level
        vsearch2seq2tax_df_selected = vsearch2seq2tax_df_selected.loc[
            vsearch2seq2tax_df_selected.rank_id >= min_tax_level_id]
        if vsearch2seq2tax_df_selected.empty:  #  no lines selected at this alignment identity threshold
            # )
            continue  #  next identity threshold
        logger.debug(
            "file: {}; line: {}; marker_name {} variant_seq {}... identity_threshold {} passed".format(__file__, inspect.currentframe().f_lineno,
                                                                       marker_name, variant_seq[1:20], identity_threshold))
        #
        # test min_tax_n
        if vsearch2seq2tax_df_selected.shape[0] < min_tax_n:
            # Logger.instance().info(
            #     "Not enought sequences are selected passing to next identity threshold."
            # )
            continue  #  next identity threshold
        logger.debug(
            "file: {}; line: {}; marker_name {} variant_seq {}... identity_threshold {} passed".format(__file__, inspect.currentframe().f_lineno,
                                                                       marker_name, variant_seq[1:20], identity_threshold))
        # continue only if selected lines
        tax_seq_id_list = vsearch2seq2tax_df_selected.tax_seq_id.tolist()
        #
        # Create lineage df
        tax_lineage_df = create_phylogenetic_line_df(tax_seq_id_list, tax_assign_sqlite)
        #
        #  Search LTG
        tax_count_perc = dataframe2ltgdefinition(tax_lineage_df)
        #
        # test min_tax_n
        if tax_count_perc.empty:
            continue  #  next identity threshold
        logger.debug(
            "file: {}; line: {}; marker_name {} variant_seq {}... identity_threshold {} passed".format(__file__, inspect.currentframe().f_lineno,
                                                                       marker_name, variant_seq[1:20], identity_threshold))
        tax_count_perc['rank_index'] = [rank_hierarchy.index(rank_name) for rank_name in
                                        tax_count_perc.index.tolist()]
        #
        #  Criteria: lineage df tax id more detailed than
        tax_count_perc.loc[tax_count_perc.rank_index >= min_tax_level_id]
        tax_count_perc_ltg = tax_count_perc.loc[tax_count_perc.rank_index >= min_tax_level_id]
        if tax_count_perc_ltg.empty:
            continue
        logger.debug(
            "file: {}; line: {}; marker_name {} variant_seq {}... identity_threshold {} passed".format(__file__, inspect.currentframe().f_lineno,
                                                                       marker_name, variant_seq[1:20], identity_threshold))
        ltg_tax_id = tax_count_perc_ltg.tax_id.tolist()[-1]
        ltg_rank_id = tax_count_perc_ltg.rank_index.tolist()[-1]
        #
        if ltg_rank_id > max_tax_resolution_id:  #  go up in lineage of ltg_tax_id up to max_tax_resolution_id
            ltg_tax_id = tax_lineage_df.loc[
                tax_lineage_df[rank_hierarchy[ltg_rank_id]] == ltg_tax_id, max_tax_resolution].unique()
        logger.debug(
            "file: {}; line: {}; marker_name {} variant_seq {}... identity_threshold {} ltg_tax_id {}".format(__file__, inspect.currentframe().f_lineno,
                                                                       marker_name, variant_seq[1:20], identity_threshold, ltg_tax_id))
        return int(ltg_tax_id)
    return int(ltg_tax_id)


def convert_fileinfo_to_otu_df(filterinfo_df):
    """
    Pivot filterinfo_df (long DF) to be used as base for the find OTU table with the tax information (wide DF)

    :param filterinfo_df: DF from marker_filterinfo_tsv in long format with variant and sample_replicate count. This are three first example lines for the MFZR marker
    variant_seq	replicate	biosample	sample_replicate	count	is_borderline	is_pseudogene_indel	is_pseudogene_codon_stop	read_average
CCTTTATCTAGTATTCGGTGCTTGGGCTGGGATAGTTGGAACAGCCCTTAGCTTACTAATCCGTGCAGAGCTTAGCCAACCTGGCGCCCTGCTCGGTGACGACCAAGTTTACAACGTGATCGTAACAGCTCATGCTTTCGTAATAATCTTCTTTATAGTAATGCCAATTATGATT	repl2	P3	P3_repl2	18	False	False	False	20
CCTTTATCTAGTATTCGGTGCTTGGGCTGGGATAGTTGGAACAGCCCTTAGCTTACTAATCCGTGCAGAGCTTAGCCAACCTGGCGCCCTGCTCGGTGACGACCAAGTTTACAACGTGATCGTAACAGCTCATGCTTTCGTAATAATCTTCTTTATAGTAATGCCAATTATGATT	repl2	P4	P4_repl2	28	False	False	False	20

    :return: Pivotted DF with one variant per line and read_count per sample_replicate in columns
    """
    # filterinfo_df.drop('Unnamed: 0', axis=1, inplace=True) # remove unnecessary column
    #
    cols = filterinfo_df.columns.tolist()
    cols.remove('biosample')
    cols.remove('replicate')
    #cols.remove('count')
    filterinfo_df = filterinfo_df.reindex(columns=cols)
    #
    #  Long to wide format for sample_replicate column
    #filterinfo_df['True'] = 1
    cols.remove('sample_replicate')
    cols.remove('count')
    filterinfo_df = pandas.pivot_table(filterinfo_df, index=cols, columns=['sample_replicate'], values='count')
    filterinfo_df.reset_index(inplace=True)
    return filterinfo_df
