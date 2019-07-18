#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import inspect
import sys

import jinja2
import yaml
import os
# import subprocess

from wopmetabarcoding.utils.constants import tempdir
from wopmetabarcoding.utils.logger import logger


def vtam_run(args_dic):
    #
    wopfile_file_name = 'Wopfile_otu.yml'
    wopfile_in_path = os.path.join(os.path.dirname(__file__), '../data', wopfile_file_name)
    wopfile_out_path = os.path.join(tempdir, wopfile_file_name)
    #
    #############################################################
    #
    # Read parameter values from param_yml file
    #
    #############################################################
    args_dic.update(yaml.load(open(args_dic['params'])))
    #
    #############################################################
    #
    # Render Wopfile
    #
    #############################################################
    with open(wopfile_in_path) as fin:
        template = jinja2.Template(fin.read())
    wopfile_rendered = template.render(args_dic)
    with open(wopfile_out_path, "w") as fout:
        fout.write(wopfile_rendered)
    args_dic['wopfile_out_path'] = wopfile_out_path
    #
    logger.debug(
        "file: {}; line: {}; Wopfile path: {}".format(__file__, inspect.currentframe().f_lineno, wopfile_out_path))
    #
    #############################################################
    #
    # Run Wopfile
    #
    #############################################################
    cmd = "wopmars -w {wopfile_out_path} -D sqlite:///{db} -v -p".format(**args_dic)
    logger.debug(
        "file: {}; line: {}; CMD: {}".format(__file__, inspect.currentframe().f_lineno, cmd))
    if args_dic['dryrun']:
        cmd = cmd + " -n"
    if args_dic['forceall']:
        cmd = cmd + " -F"
    if 'targetrule' is args_dic:
        cmd = cmd + " -t {targetrule}".format(**args_dic)
    # import pdb; pdb.set_trace()
    # p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    os.system(cmd)
    # stdout = p.communicate()
    # p.wait()  # wait program to finish
    sys.exit(0)

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--db', dest='db', nargs=1, help="SQLITE file with DB", required=True)
    parser.add_argument('--fastainfo', dest='fastainfo', nargs=1, help="TSV file with FASTA sample information", required=True)
    parser.add_argument('--fastadir', dest='fastadir', nargs=1, help="Directory with FASTA files", required=True)
    parser.add_argument('--outdir', nargs=1, help="Directory for output", required=True)
    parser.add_argument('--params', nargs=1, help="YML file with parameter values", required=True)
    parser.add_argument('--filter_lfn_variant', nargs=1, help="Boolean 0|1 to filter_lfn_variant (1) or filter_lfn_variant_replicate (0)", required=True)
    parser.add_argument('-F', '--forceall', action='store_true', help="Force argument of WopMars", required=False)
    parser.add_argument('-t', '--targetrule', nargs=1, help="Execute the workflow to the given RULE: SampleInformation, ...", required=False)
    parser.add_argument('-n', '--dry-run', dest='dryrun', action='store_true', help="Only display what would have been done.")
    return parser

def main():
    parser = create_parser()
    args = parser.parse_args()
    #
    args_dic = {
        'db': args.db[0],
        'fastainfo': args.fastainfo[0],
        'fastadir': args.fastadir[0],
        'outdir': args.outdir[0],
        'filter_lfn_variant': args.filter_lfn_variant[0],
        'params': args.params[0],
        'forceall': args.forceall,
        'sortreads': os.path.join(args.outdir[0], 'sortreads.tsv'),
        'outtable': os.path.join(args.outdir[0], 'otutable.tsv'),
        'dryrun': args.dryrun,
    }
    if not args.targetrule is None:
        args_dic['targetrule'] = args.targetrule[0]
    vtam_run(args_dic)

if __name__=='__main__':
    main()


