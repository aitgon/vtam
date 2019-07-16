#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

import jinja2
import os
import subprocess

from wopmetabarcoding.utils.constants import tempdir

def vtam_run(args_dic):
    #
    wopfile_in_path = os.path.join(os.path.dirname(__file__), '../data', 'Wopfile_merge.yml')
    wopfile_out_path = os.path.join(tempdir, 'Wopfile_merge.yml')
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
    print(wopfile_out_path)
    #
    #############################################################
    #
    # Run Wopfile
    #
    #############################################################
    # cmd = "wopmars -w {} -D sqlite:///{} -v -p".format(wopfile_out_path, wopdb)
    cmd = "wopmars -w {wopfile_out_path} -D sqlite:///{wopdb} -v -p".format(**args_dic)
    if args_dic['forceall']:
        cmd = cmd + " -F"
    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    stdout = p.communicate()
    p.wait()  # wait program to finish
    sys.exit(0)

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--wopdb', dest='wopdb', nargs=1, help="SQLITE file with WopMars DB")
    parser.add_argument('--fastqinfo', dest='fastqinfo', nargs=1, help="TSV file with FASTQ sample information")
    parser.add_argument('--fastainfo', dest='fastainfo', nargs=1, help="TSV file with FASTA sample information")
    parser.add_argument('--fastqdir', dest='fastqdir', nargs=1, help="Directory with FASTQ files")
    parser.add_argument('--fastadir', dest='fastadir', nargs=1, help="Directory with FASTA files")
    parser.add_argument('-F', '--forceall', action='store_true', help="Force argument of WopMars")
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
        'forceall': args.forceall
    }
    # vtam_run(args.wopdb[0], args.fastqinfo[0], args.fastainfo[0], args.fastqdir[0], args.fastadir[0], args.forceall)
    vtam_run(args_dic)

if __name__=='__main__':
    main()

