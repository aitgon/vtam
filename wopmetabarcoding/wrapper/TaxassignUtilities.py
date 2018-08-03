import inspect

import pandas, os, sqlite3, itertools, csv
from Bio import SeqIO
from numpy import nan
from wopmars.utils.Logger import Logger

from wopmetabarcoding.utils.constants import tempdir
from wopmetabarcoding.utils.logger import logger, LOGGER_LEVEL
from wopmetabarcoding.utils.PathFinder import PathFinder

rank_hierarchy =['no rank', 'phylum', 'superclass', 'class', 'subclass', 'infraclass', 'superorder', 'order', 'suborder', 'infraorder', 'family', 'subfamily', 'genus', 'subgenus', 'species', 'subspecies']

def f_taxid2taxname(tax_id_list, tax_assign_sqlite):
    con = sqlite3.connect(tax_assign_sqlite)
    taxid2taxname_dic = {}
    for tax_id in tax_id_list:
        tax_name = nan
        sql = "SELECT tax_name FROM seq2tax2parent WHERE tax_id = ?"
        cur = con.cursor()
        cur.execute(sql, (tax_id,))
        row = cur.fetchone()
        if not row is None:
            tax_name = str(row[0])
        cur.close()
        taxid2taxname_dic[tax_id] = tax_name
    con.close()
    return taxid2taxname_dic

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
    Given a list of taxon sequence ids (tax_seq_id_list), create a df with a taxon lineage per df

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
        try:
            tax_id = int(row[0])
        except:
            tax_id = row[0]
        rank_name = row[2]
        tax_parent_id = int(row[3])
        # while not row is None or row[0] != 1:
        while row[0] != 1:
            tax_lineage[rank_name] = tax_id
            filter_string = tax_parent_id
            sql = "SELECT tax_id, tax_name, parent_id, rank FROM tax2parent WHERE tax_id = ?"
            cur.execute(sql, (filter_string,))
            row = cur.fetchone()
            if row is None:
                break
            tax_id = int(row[0])
            rank_name = row[3]
            tax_parent_id = int(row[2])
        lineage_list.append(tax_lineage)
        cur.close()
    conn.close()
    tax_lineage_df = pandas.DataFrame(lineage_list)
    tax_lineage_df = tax_lineage_df[tax_lineage_header]
    #
    tax_lineage_df.fillna(value=nan, inplace=True)
    tax_lineage_df.replace(to_replace='', value=nan, inplace=True)
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


def f_majoritytaxid2percentage(tax_lineage_df):
    """
    Given the lineages for a list of tax_seq_ids, this function will compute the percentage of the most common tax_id at each rank level

    :param tax_lineage_df: Output of tax_lineage_df that looks like this
  tax_seq_id  no rank  phylum  class  subclass  infraclass  order  suborder  \
6      3433568   131567    6656  50557      7496       33340   7524       NaN
9      4000774   131567    6656  50557      7496       33340   7088       NaN
17     4472326   131567    6656  50557      7496       33340   7088   41191.0
7      4472342   131567    6656  50557      7496       33340   7088   41191.0
13     4473884   131567    6656  50557      7496       33340   7088   41191.0

    :return:a df with this information

(Pdb) tax_count_perc
           tax_id  count   perc
phylum     6656.0     20  100.0
class     50557.0     20  100.0
subclass   7496.0     20  100.0
    """
    tax_count_perc = pandas.DataFrame({'tax_id': tax_lineage_df.apply(lambda x: x.value_counts().index[0], axis=0)})
    tax_count_perc['count'] = tax_lineage_df.apply(lambda x: x.value_counts().iloc[0], axis=0)
    tax_count_perc['perc'] = tax_count_perc['count'] / tax_lineage_df.shape[0] * 100
    tax_count_perc.drop(['tax_seq_id', 'no rank'], axis=0, inplace=True)
    # converting to int, makes it sure that there are non empty cells
    tax_count_perc[['tax_id', 'count']] = tax_count_perc[['tax_id', 'count']].astype('int')
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
    sql = "SELECT target, identity_thresold FROM alignedtsv WHERE query_variant == '{}'".format(variant_seq)
    vsearch_output_for_variant_df = pandas.read_sql(sql=sql, con=con)
    con.close()
    return vsearch_output_for_variant_df


def f_taxlineage_to_ltg(tax_lineage_df, max_tax_resolution_id):
    #
    # 2. Select majority taxon_id at each level
    # 3. Compute percentage of majority taxon_id at each level
    tax_count_perc = f_majoritytaxid2percentage(tax_lineage_df)
    #
    if tax_count_perc.empty:
        return None  # next identity threshold
    # 4. Select majority taxon_id with more 90% presence at a given level
    tax_count_perc = tax_count_perc.loc[tax_count_perc.perc >= 90.0]
    #
    if tax_count_perc.empty:
        return None  # next identity threshold
    tax_count_perc['rank_index'] = [rank_hierarchy.index(rank_name) for rank_name in tax_count_perc.index.tolist()]
    #
    # 5. The rank level of the ltg must be less detailed than max_tax_resolution
    tax_count_perc = tax_count_perc.loc[tax_count_perc.rank_index <= max_tax_resolution_id]
    #
    # This line was commented out but be careful, because unclear if really unnecessary.
    # tax_count_perc_ltg = tax_count_perc.loc[tax_count_perc.rank_index >= min_tax_level_id]
    if tax_count_perc.empty:
        return None
    #
    # 6. The LTG is the most detailed tax_id among the remaining tax_id
    ltg_tax_id = tax_count_perc.tail(1)['tax_id'].values[0]
    return ltg_tax_id

def f_variant_vsearch_output_to_ltg(variant_id, vsearch_output_for_variant_df_pkl, tax_assign_sqlite, tax_assign_pars_tsv):
    ltg_tax_id = nan # default ltg_tax_id
    vsearch_output_for_variant_df = pandas.read_pickle(vsearch_output_for_variant_df_pkl, compression='gzip')
    vsearch_output_for_variant_df.columns = ["tax_seq_id", "alignment_identity"]
    #
    tax_seq_id_list = vsearch_output_for_variant_df.tax_seq_id.tolist()
    #
    seq2tax_df = seq2tax_db_sqlite_to_df(tax_assign_sqlite, tax_seq_id_list)
    #
    # tax_assign_pars df
    #names = ["identity_threshold", "min_tax_level", "max_tax_resolution", "min_tax_n"]
    # header is now provided in the pars tsv file
    tax_assign_pars_df = pandas.read_csv(tax_assign_pars_tsv, sep="\t", header=0)
    #
    # Merge of the vsearch alignment, the sequence and taxa information
    vsearch_output_for_variant_df[["tax_seq_id"]] = vsearch_output_for_variant_df[["tax_seq_id"]].astype('int64')
    seq2tax_df[["tax_seq_id"]] = seq2tax_df[["tax_seq_id"]].astype('int64')
    vsearch_output_for_variant_df = pandas.merge(vsearch_output_for_variant_df, seq2tax_df, left_on="tax_seq_id", right_on="tax_seq_id")
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
        min_tax_level_id = rank_hierarchy.index(min_tax_level) # id of tax_level
        max_tax_resolution_id = rank_hierarchy.index(max_tax_resolution)
        #
        ########################
        # Aucun hit a 100 => on passe au niveau de similarite suivant
        # At given threshold, assignation fails if no hits
        ########################
        vsearch2seq2tax_df_selected = vsearch_output_for_variant_df.loc[
            vsearch_output_for_variant_df.alignment_identity >= identity_threshold]
        logger.debug(
            "file: {}; line: {}; variant_id {}... identity_threshold {} vsearch2seq2tax_df_selected shape {}".format(__file__, inspect.currentframe().f_lineno,
                                                                                                                                     variant_id, identity_threshold, vsearch2seq2tax_df_selected.shape))
        if vsearch2seq2tax_df_selected.empty:  # no lines selected at this alignment identity threshold
            continue  # next identity threshold
        logger.debug(
            "file: {}; line: {}; variant_id {}... identity_threshold {} passed".format(__file__, inspect.currentframe().f_lineno,
                                                                                                       variant_id, identity_threshold))
        #
        ########################################################################
        # A ce pourcentage on accepte que des hits qui sont annotes au niveau famille (min_tax_level_id) ou plus precise (rank_id plus eleve) (genus, espece) => On garde les 6 hits
        # At given threshold, discard hits that are annotated with tax_id less detailed than min_tax_level_id
        # Assignation fails if no hits after selection
        ########################################################################
        vsearch2seq2tax_df_selected = vsearch2seq2tax_df_selected.loc[vsearch2seq2tax_df_selected.rank_id >= min_tax_level_id]
        if vsearch2seq2tax_df_selected.empty:  # no lines selected at this alignment identity threshold
            # )
            continue  # next identity threshold
        logger.debug(
            "file: {}; line: {}; variant_id {}... identity_threshold {} passed".format(__file__, inspect.currentframe().f_lineno,
                                                                                                       variant_id, identity_threshold))
        #
        ########################################################################
        # Il faut au moins 3 (min_tax_n) taxa parmi les hits => on a 3 (86610, 6115, 6116), donc c'est bon (si non, on passe au niveau de pourcentage superieure)
        # At given threshold, assignation fails if the number of different taxa is strictly lower than min_tax_n
        ########################################################################
        if vsearch2seq2tax_df_selected.shape[0] < min_tax_n:
            continue  # next identity threshold
        logger.debug(
            "file: {}; line: {}; variant_id {}... identity_threshold {} passed".format(__file__, inspect.currentframe().f_lineno,
                                                                                                       variant_id, identity_threshold))
        # continue only if selected lines
        tax_seq_id_list = vsearch2seq2tax_df_selected.tax_seq_id.tolist()
        #
        ########################################################################
        # On prend le plus petit groupe taxonomique qui contient au moins 90% des hits
        # We have to
        # 1. Create tax_lineage_df
        # 2. Select majority taxon_id at each level
        # 3. Compute percentage of majority taxon_id at each level
        # 4. Select majority taxon_id with more 90% presence at a given level
        # 5. Select tax_id's with rank level be less detailed than max_tax_resolution
        # 6. The LTG is the most detailed tax_id among the remaining tax_id
        ########################################################################
        #
        # Create lineage df
        logger.debug(
            "file: {}; line: {}; create_phylogenetic_line_df".format(__file__, inspect.currentframe().f_lineno))
        # 1. Create tax_lineage_df
        tax_lineage_df = create_phylogenetic_line_df(tax_seq_id_list, tax_assign_sqlite)
        # 2. Select majority taxon_id at each level
        # 3. Compute percentage of majority taxon_id at each level
        # 4. Select majority taxon_id with more 90% presence at a given level
        # 5. Select tax_id's with rank level be less detailed than max_tax_resolution
        # 6. The LTG is the most detailed tax_id among the remaining tax_id
        ltg_tax_id = f_taxlineage_to_ltg(tax_lineage_df, max_tax_resolution_id)
        if ltg_tax_id is None:
            continue  # next identity threshold
        else:
            return ltg_tax_id
    return nan

def f_variant_vsearch_output_to_ltg_bak(variant_seq, marker_name, vsearch_output_for_variant_df_pkl, tax_assign_sqlite, tax_assign_pars_tsv):
    ltg_tax_id = nan # default ltg_tax_id
    vsearch_output_for_variant_df = pandas.read_pickle(vsearch_output_for_variant_df_pkl, compression='gzip')
    vsearch_output_for_variant_df.columns = ["tax_seq_id", "alignment_identity"]
    #
    tax_seq_id_list = vsearch_output_for_variant_df.tax_seq_id.tolist()
    #
    seq2tax_df = seq2tax_db_sqlite_to_df(tax_assign_sqlite, tax_seq_id_list)
    #
    # tax_assign_pars df
    #names = ["identity_threshold", "min_tax_level", "max_tax_resolution", "min_tax_n"]
    # header is now provided in the pars tsv file
    tax_assign_pars_df = pandas.read_csv(tax_assign_pars_tsv, sep="\t", header=0)
    #
    # Merge of the vsearch alignment, the sequence and taxa information
    vsearch_output_for_variant_df[["tax_seq_id"]] = vsearch_output_for_variant_df[["tax_seq_id"]].astype('int64')
    seq2tax_df[["tax_seq_id"]] = seq2tax_df[["tax_seq_id"]].astype('int64')
    vsearch_output_for_variant_df = pandas.merge(vsearch_output_for_variant_df, seq2tax_df, left_on="tax_seq_id", right_on="tax_seq_id")
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
        min_tax_level_id = rank_hierarchy.index(min_tax_level) # id of tax_level
        max_tax_resolution_id = rank_hierarchy.index(max_tax_resolution)
        #
        ########################
        # Aucun hit a 100 => on passe au niveau de similarite suivant
        # At given threshold, assignation fails if no hits
        ########################
        vsearch2seq2tax_df_selected = vsearch_output_for_variant_df.loc[
            vsearch_output_for_variant_df.alignment_identity >= identity_threshold]
        logger.debug(
            "file: {}; line: {}; marker_name {} variant_seq {}... identity_threshold {} vsearch2seq2tax_df_selected shape {}".format(__file__, inspect.currentframe().f_lineno,
            marker_name, variant_seq[1:20], identity_threshold, vsearch2seq2tax_df_selected.shape))
        if vsearch2seq2tax_df_selected.empty:  # no lines selected at this alignment identity threshold
            continue  # next identity threshold
        logger.debug(
            "file: {}; line: {}; marker_name {} variant_seq {}... identity_threshold {} passed".format(__file__, inspect.currentframe().f_lineno,
                                                                       marker_name, variant_seq[1:20], identity_threshold))
        #
        ########################################################################
        # A ce pourcentage on accepte que des hits qui sont annotes au niveau famille (min_tax_level_id) ou plus precise (rank_id plus eleve) (genus, espece) => On garde les 6 hits
        # At given threshold, discard hits that are annotated with tax_id less detailed than min_tax_level_id
        # Assignation fails if no hits after selection
        ########################################################################
        vsearch2seq2tax_df_selected = vsearch2seq2tax_df_selected.loc[vsearch2seq2tax_df_selected.rank_id >= min_tax_level_id]
        if vsearch2seq2tax_df_selected.empty:  # no lines selected at this alignment identity threshold
            # )
            continue  # next identity threshold
        logger.debug(
            "file: {}; line: {}; marker_name {} variant_seq {}... identity_threshold {} passed".format(__file__, inspect.currentframe().f_lineno,
                                                                       marker_name, variant_seq[1:20], identity_threshold))
        #
        ########################################################################
        # Il faut au moins 3 (min_tax_n) taxa parmi les hits => on a 3 (86610, 6115, 6116), donc c'est bon (si non, on passe au niveau de pourcentage superieure)
        # At given threshold, assignation fails if the number of different taxa is strictly lower than min_tax_n
        ########################################################################
        if vsearch2seq2tax_df_selected.shape[0] < min_tax_n:
            continue  # next identity threshold
        logger.debug(
            "file: {}; line: {}; marker_name {} variant_seq {}... identity_threshold {} passed".format(__file__, inspect.currentframe().f_lineno,
                                                                       marker_name, variant_seq[1:20], identity_threshold))
        # continue only if selected lines
        tax_seq_id_list = vsearch2seq2tax_df_selected.tax_seq_id.tolist()
        #
        ########################################################################
        # On prend le plus petit groupe taxonomique qui contient au moins 90% des hits
        # We have to
        # 1. Create tax_lineage_df
        # 2. Select majority taxon_id at each level
        # 3. Compute percentage of majority taxon_id at each level
        # 4. Select majority taxon_id with more 90% presence at a given level
        # 5. Select tax_id's with rank level be less detailed than max_tax_resolution
        # 6. The LTG is the most detailed tax_id among the remaining tax_id
        ########################################################################
        #
        # Create lineage df
        logger.debug(
            "file: {}; line: {}; create_phylogenetic_line_df".format(__file__, inspect.currentframe().f_lineno))
        # 1. Create tax_lineage_df
        tax_lineage_df = create_phylogenetic_line_df(tax_seq_id_list, tax_assign_sqlite)
        # 2. Select majority taxon_id at each level
        # 3. Compute percentage of majority taxon_id at each level
        # 4. Select majority taxon_id with more 90% presence at a given level
        # 5. Select tax_id's with rank level be less detailed than max_tax_resolution
        # 6. The LTG is the most detailed tax_id among the remaining tax_id
        ltg_tax_id = f_taxlineage_to_ltg(tax_lineage_df, max_tax_resolution_id)
        if ltg_tax_id is None:
            continue  # next identity threshold
        else:
            return ltg_tax_id
    return nan


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
    filterinfo_df = pandas.pivot_table(filterinfo_df, index=cols, columns=['sample_replicate'], values='count', fill_value=0)
    filterinfo_df.reset_index(inplace=True)
    return filterinfo_df
