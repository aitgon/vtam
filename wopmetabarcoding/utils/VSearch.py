# 'vsearch –uchime_denovo ' + sorted_repl_fasta_name + ' --borderline ' + borderline_repl_filename +
#     ' --nonchimeras '+ nonchimera_repl_filename + ' --chimeras ' + chimera_repl_filename, shell=True


import os
import subprocess

class VSearch:

    def __init__(self):
        self.params = {
        }

    def run(self):
        cmd_args = [str(item) for item in list(sum(list(self.params.items()), ()))]
        cmd_line = ['vsearch'] + cmd_args
        with open(os.devnull, 'w') as FNULL:
            p = subprocess.Popen(cmd_line, stdout=FNULL, stderr=FNULL)
            p.wait()

class VSearch1(VSearch):

    def __init__(self, db, usearch_global, id=None, maxhits=None, maxrejects=None, maxaccepts=None, minseqlength=None, userout=None, userfields=None):
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


def main():
    vsearch = VSearch()
    vsearch.run()
    vsearch_params = {'db' : 'dbpath', 'usearch_global' : 'usearch_global', 'maxrejects': 1}
    vsearch1 = VSearch1(**vsearch_params)
    vsearch1.run()

if __name__ == '__main__':
    main()


