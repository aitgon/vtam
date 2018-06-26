import pandas, os, sqlite3, itertools, csv
from Bio import SeqIO
from numpy import nan
from wopmetabarcoding.utils.constants import tempdir

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
    os.system(
        "vsearch --makeudb_usearch " + taxassign_db_fasta + " --output " + udb_database
    )


def alignment_vsearch(filtered_variants_fasta, taxassign_db_fasta, output_tsv):
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


def sub_fasta_creator(marker_variant_fasta, sequence_number, marker_name):
    """
    Split a FASTA file into FASTA pieces of given sequence number

    :param marker_variant_fasta: Path to FASTA file with marker variants
    :param sequence_number: Maximal sequence number in split FASTA files
    :param marker_name: Marker name to prefix split FASTA files
    :return: List with paths to split FASTA files
    """
    record_iter = SeqIO.parse(open(marker_variant_fasta), "fasta")
    print(str(record_iter))
    sub_fasta_path_list = []
    for i, batch in enumerate(batch_iterator(record_iter, sequence_number)):
        sub_fasta_path = "group_{}_{}.fasta".format(i + 1, marker_name)
        # filename = filename + "_" + marker_name + ".fasta"
        sub_fasta_path = os.path.join(tempdir, sub_fasta_path)
        with open(sub_fasta_path, "w") as handle:
            count = SeqIO.write(batch, handle, "fasta")
        sub_fasta_path_list.append(sub_fasta_path)
    return sub_fasta_path_list


def create_tsv_per_variant(filename, db_to_create):
    conn = sqlite3.connect(db_to_create)
    conn.execute("DROP TABLE IF EXISTS alignedtsv")
    cur = conn.cursor()
    conn.execute(
        "CREATE TABLE alignedtsv (id INTEGER PRIMARY KEY AUTOINCREMENT, query_variant VARCHAR, target VARCHAR, identity_thresold FLOAT)"
    )
    with open(filename, 'r') as fin:
        print(filename)
        for line in fin:
            line = line.strip().split("\t")
            cur.execute("INSERT INTO alignedtsv (query_variant, target, identity_thresold) VALUES (?,?,?);", line)
    cur.close()
    conn.commit()
    conn.close()




