#!/usr/bin/env python
# -*- coding: utf-8 -*-

from vtam.utils.TaxonomyDB import TaxonomyDB

import argparse

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', dest='output', action='store', help="Path to sqlite DB with taxonomy information",
                        required=True)
    # parser.add_argument('--datadir', action='store',
    #         help="Directory that will be used for __download_ncbi_taxonomy_dump NCBI data. If this argument is absent, "
    #              "the data will be downloaded to the temporary directory.", required=False)
    parser.add_argument('--precomputed', dest='precomputed', action='store_true', default=True,
                             help="Should we use precomputed and old taxonomy DB or new and potentially unstable one", required=False)
    return parser

def main():
    parser = create_parser()
    args = parser.parse_args()
    taxonomydb = TaxonomyDB(output=vars(args)['output'], precomputed=vars(args)['precomputed'], )
    taxonomydb.create_taxonomy_db()


if __name__=='__main__':
    main()
