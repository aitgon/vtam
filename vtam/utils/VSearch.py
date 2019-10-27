import inspect
import os
import subprocess

from vtam.utils.Logger import Logger

class VSearch2():

    def __init__(self, parameters):
        """Creates a vsearch object based on the parameters

        :param parameters: Dictionary with key:value for store parameters or key:None for other parameters
        :return: void
        """
        self.parameters = parameters

    def run(self):
        command = 'vsearch'
        for param in self.parameters:
            if not self.parameters[param] is None:
                command += ' {} {}'.format(param, self.parameters[param])
            else:
                command += ' {}'.format(param)
        Logger.instance().info(command)
        run_result = subprocess.run(command.split(), stdout=subprocess.PIPE)
        Logger.instance().info(run_result.stdout)


class VSearch:

    def __init__(self):
        self.params = {
        }

    def run(self):
        cmd_args = [str(item) for item in list(sum(list(self.params.items()), ()))]
        # TODO: include threads argument
        cmd_line = ['vsearch'] + cmd_args
        with open(os.devnull, 'w') as FNULL:
            Logger.instance().debug(
                "file: {}; line: {}; Vsearch command: {}".format(__file__, inspect.currentframe().f_lineno, " ".join(cmd_line)))
            p = subprocess.Popen(cmd_line, stdout=FNULL, stderr=FNULL)
            p.wait()


class VSearchUsearchGlobal(VSearch):

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



class VsearchChimera(VSearch):

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
    vsearch1 = VSearchUsearchGlobal(**vsearch_params)
    vsearch1.run()


if __name__ == '__main__':
    main()


