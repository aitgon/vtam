#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import jinja2
import os
import subprocess

from wopmetabarcoding.utils.constants import tempdir

def vtam_run(wopdb, fastainfo, fastadir):
    #
    wopfile_in_path = os.path.join(os.path.dirname(__file__), '../data', 'Wopfile.yml')
    wopfile_out_path = os.path.join(tempdir, 'Wopfile.yml')
    #
    #############################################################
    #
    # Render Wopfile
    #
    #############################################################
    with open(wopfile_in_path) as fin:
        template = jinja2.Template(fin.read())
    context = {'fastainfo' : fastainfo, 'fastadir' : fastadir}
    wopfile_rendered = template.render(context)
    with open(wopfile_out_path, "w") as fout:
        fout.write(wopfile_rendered)
    #
    #############################################################
    #
    # Run Wopfile
    #
    #############################################################
    cmd = "wopmars -w {} -D sqlite:///{} -v -p".format(wopfile_out_path, wopdb)
    # import pdb; pdb.set_trace()
    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    stdout = p.communicate()
    p.wait()  # wait program to finish

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--wopdb', dest='wopdb', nargs=1, help="SQLITE file with WopMars DB")
    parser.add_argument('--fastainfo', dest='fastainfo', nargs=1, help="TSV file with FASTA sample information")
    parser.add_argument('--fastadir', dest='fastadir', nargs=1, help="Directory with FASTA files")
    return parser

def main():
    parser = create_parser()
    args = parser.parse_args()
    #
    vtam_run(args.wopdb[0], args.fastainfo[0], args.fastadir[0])

if __name__=='__main__':
    main()

