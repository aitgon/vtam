#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

import jinja2
import os
import subprocess

import yaml

from pathlib import Path

from wopmetabarcoding.utils.OptionManager import OptionManager
from wopmetabarcoding.utils.utilities import tempdir
from wopmetabarcoding.utils.PathManager import PathFinder


def vtam_run(args_dic):
    #
    #
    #############################################################
    #
    # Read parameter values from param_yml file
    #
    #############################################################
    args_dic.update(yaml.load(open(args_dic['params']), Loader=yaml.SafeLoader))    #
    wopfile_file_name = 'Wopfile_merge.yml'
    wopfile_in_path = os.path.join(os.path.dirname(__file__), '../data', wopfile_file_name)
    wopfile_out_path = os.path.join(tempdir, wopfile_file_name)
    #
    #############################################################
    #
    # Render Wopfile
    #
    #############################################################
    with open(wopfile_in_path) as fin:
        template = jinja2.Template(fin.read())
    # context = {'fastqinfo' : fastqinfo, 'fastainfo' : fastainfo, 'fastqdir' : fastqdir, 'fastadir' : fastadir}
    wopfile_rendered = template.render(args_dic)
    with open(wopfile_out_path, "w") as fout:
        fout.write(wopfile_rendered)
    args_dic['wopfile_out_path'] = wopfile_out_path
    #
    #############################################################
    #
    # Run Wopfile
    #
    #############################################################
    # cmd = "wopmars -w {} -D sqlite:///{} -v -p".format(wopfile_out_path, wopdb)
    cmd = "wopmars -w {wopfile_out_path} -D sqlite:///{wopdb} -p -v".format(**args_dic)
    if args_dic['dryrun']:
        cmd = cmd + " -n"
    if args_dic['forceall']:
        cmd = cmd + " -F"
    if not args_dic['log'] is None:
        PathFinder.mkdir_p(os.path.dirname(args_dic['log']))
        Path(args_dic['log']).touch() # touch log
    #     cmd = cmd + " --log " + args_dic['log']
    from wopmetabarcoding.utils.logger import logger
    logger.info(cmd)
    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    stdout = p.communicate()
    p.wait()  # wait program to finish
    sys.exit(0)

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--wopdb', dest='wopdb', nargs=1, help="SQLITE file with WopMars DB", required=True)
    parser.add_argument('--fastqinfo', dest='fastqinfo', nargs=1, help="TSV file with FASTQ sample information",
                        required=True, type=os.path.abspath)
    parser.add_argument('--fastainfo', dest='fastainfo', nargs=1, help="TSV file with FASTA sample information",
                        required=True, type=os.path.abspath)
    parser.add_argument('--fastqdir', dest='fastqdir', nargs=1, help="Directory with FASTQ files", required=True,
                        type=os.path.abspath)
    parser.add_argument('--fastadir', dest='fastadir', nargs=1, help="Directory with FASTA files", required=True,
                        type=os.path.abspath)
    parser.add_argument('--params', nargs=1, help="YML file with parameter values", required=True)
    parser.add_argument('-F', '--forceall', action='store_true', help="Force argument of WopMars", required=False)
    parser.add_argument('-n', '--dry-run', dest='dryrun', action='store_true',
                        help="Only display what would have been done.")
    parser.add_argument('--log', nargs=1, dest='log',
                        help="Write logs in FILE file [default: $HOME/.wopmars/wopmars.log].", required=False,
                        type=os.path.abspath)
    return parser

def main():
    parser = create_parser()
    args = parser.parse_args()
    #
    args_dic = {
        'wopdb': args.wopdb[0],
        'fastqinfo': args.fastqinfo[0],
        'fastainfo': args.fastainfo[0],
        'fastqdir': args.fastqdir[0],
        'fastadir': args.fastadir[0],
        'forceall': args.forceall,
        'dryrun': args.dryrun,
        'params': args.params[0],
        'log': None,
    }
    if not args.log is None:
        OptionManager.instance()['--log'] = args.log[0]
        args_dic['log'] = args.log[0]
    vtam_run(args_dic)

if __name__=='__main__':
    main()

