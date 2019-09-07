import inspect
import os
import subprocess

import sys

from wopmetabarcoding.utils.logger import logger


class VSearch:

    def __init__(self):
        self.params = {
        }

    def run(self):
        cmd_args = [str(item) for item in list(sum(list(self.params.items()), ()))]
        cmd_line = ['vsearch'] + cmd_args
        with open(os.devnull, 'w') as FNULL:
            logger.debug(
                "file: {}; line: {}; Vsearch command: {}".format(__file__, inspect.currentframe().f_lineno, " ".join(cmd_line)))
            p = subprocess.Popen(cmd_line, stdout=FNULL, stderr=FNULL)
            p.wait()


class VSearch1(VSearch):

    def __init__(self, db, usearch_global, id=None, maxhits=None, maxrejects=None, maxaccepts=None, minseqlength=None, userout=None, userfields=None, query_cov=None):
        #
        VSearch.__init__(self)
        self.params["--db"] = db
        self.params["--usearch_global"] = usearch_global
        #
        if not id is None: self.params['--id'] = id
        if not maxhits is None: self.params['--maxhits'] = maxhits
        if not maxrejects is None: self.params['--maxrejects'] = maxrejects
        if not maxaccepts is None: self.params['--maxaccepts'] = maxaccepts
        if not minseqlength is None: self.params['--minseqlength'] = minseqlength
        if not userout is None: self.params['--userout'] = userout
        if not userfields is None: self.params['--userfields'] = userfields
        if not query_cov is None: self.params['--query_cov'] = query_cov


class Vsearch2(VSearch):

    def __init__(self, sortbysize, output):
        VSearch.__init__(self)
        self.params["--sortbysize"] = sortbysize
        self.params["--output"] = output


class Vsearch3(VSearch):

    def __init__(self, uchime_denovo, borderline, nonchimeras, chimeras):
        VSearch.__init__(self)
        self.params["--uchime_denovo"] = uchime_denovo
        self.params["--borderline"] = borderline
        self.params["--nonchimeras"] = nonchimeras
        self.params["--chimeras"] = chimeras


def main():
    vsearch = VSearch()
    vsearch.run()
    vsearch_params = {'db' : 'dbpath', 'usearch_global' : 'usearch_global', 'maxrejects': 1}
    vsearch1 = VSearch1(**vsearch_params)
    vsearch1.run()


if __name__ == '__main__':
    main()


