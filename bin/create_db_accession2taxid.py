#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import gzip
import inspect
import shutil
import sqlite3
import urllib.request


import errno, numpy, os, pandas, urllib

import os

from wopmetabarcoding.utils.constants import VTAM_DATA_DIR
from wopmetabarcoding.utils.logger import logger



def download():
    try:
        os.makedirs(VTAM_DATA_DIR)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    # Download files
    remotefile = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz"
    nucl_gb_accession2taxid_gz = os.path.join(VTAM_DATA_DIR, os.path.basename(remotefile))
    logger.debug(
        "file: {}; line: {}; Downloading nucl_gb_accession2taxid_gz".format(__file__, inspect.currentframe().f_lineno))
    if not os.path.isfile(nucl_gb_accession2taxid_gz): # TODO verify MD5
        urllib.request.urlretrieve(remotefile, nucl_gb_accession2taxid_gz)
    return nucl_gb_accession2taxid_gz

def f_create_nucl_gb_accession2taxid_db(accession2tax_id_sqlite_path):

    new_taxdump_path = download()

    nucl_gb_accession2taxid_tsv = os.path.join(VTAM_DATA_DIR, "nucl_gb.accession2taxid.tsv")
    if not (os.path.isfile(nucl_gb_accession2taxid_tsv)):
        logger.debug(
            "file: {}; line: {}; Gunziping nucl_gb.accession2taxid.gz".format(__file__,inspect.currentframe().f_lineno))
        with gzip.open(new_taxdump_path, 'rb') as f_in:
            with open(nucl_gb_accession2taxid_tsv, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    #
    # Import nucl_gb.accession2taxid.gz to sqlite

    # nucl_gb_accession2taxid_sqlite = nucl_gb_accession2taxid_tsv.replace("tsv", "sqlite")
    logger.debug(
        "file: {}; line: {}; Import nucl_gb.accession2taxid.gz to sqlite".format(__file__,
                                                                                 inspect.currentframe().f_lineno))
    conn = sqlite3.connect(accession2tax_id_sqlite_path)
    conn.execute("DROP TABLE IF EXISTS nucl_gb_accession2taxid")
    cur = conn.cursor()
    conn.execute(
        "CREATE TABLE nucl_gb_accession2taxid (gb_accession VARCHAR, tax_id INTEGER, unique(gb_accession, tax_id))"
    )
    cur.close()
    conn.commit()
    chunksize = 10 ** 6
    for chunk_df in pandas.read_csv(os.path.join(VTAM_DATA_DIR, "nucl_gb.accession2taxid.tsv"), sep="\t",
                                    usecols=[1, 2],
                                    names=['gb_accession', 'tax_id'], chunksize=chunksize, header=0):
        chunk_df.to_sql(name="nucl_gb_accession2taxid", con=conn, if_exists='append', index=False)
    conn.close()

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', dest='output', help="Path to sqlite DB with accession2tax_id information")
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()
    f_create_nucl_gb_accession2taxid_db(args.output)

if __name__=='__main__':
    main()
