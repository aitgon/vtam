import tempfile

import jinja2
import os

import sys

from wopmetabarcoding.utils.Singleton import Singleton
from wopmetabarcoding.utils.constants import parameters_numerical
from wopmetabarcoding.utils.utilities import tempdir


class WopmarsRunner(Singleton):

    def __init__(self, subcommand, parameters):
        """

        :param subcommand: takes one of three values: merge, otu or optimize
        :param parameters: dictionnary (OptionManager.instance()) with command
        """
        self.subcommand = subcommand
        ##################################
        #
        # Load default numerical parameters and overwrite with custom parameters
        #
        ##################################
        self.parameters = parameters_numerical.copy()
        for k in parameters:
            self.parameters[k] = parameters[k]
        #
        self.wopfile_path = None


    def create_wopfile(self, path=None):
        """

        :param wopfile_path: Path of output wopfile
        :return: tuple (wopfile_path, wopfile_content)
        """
        # Create path
        wopfile_path = path
        if wopfile_path is None:
            wopfile_path = tempfile.NamedTemporaryFile().name
        self.wopfile_path = wopfile_path
        #
        # wopfile_file_name = 'Wopfile_{}.yml'.format(self.subcommand)
        wopfile_template_path\
            = os.path.join(os.path.dirname(__file__), '../../data/Wopfile_{subcommand}.yml'.format(subcommand=self.subcommand))
        with open(wopfile_template_path) as fin:
            template = jinja2.Template(fin.read())
        wopfile_content = template.render(self.parameters)
        ################
        #
        # return
        #
        ################
        with open(wopfile_path, "w") as fout:
            fout.write(wopfile_content)
        return (wopfile_path, wopfile_content)


    def get_wopmars_command(self):
        """

        :param wopfile_out_path: Path of output wopfile
        :return: string with output path of wopfile
        """

        ###################
        #
        # Base wopmars command
        #
        ###################
        if self.wopfile_path is None:
            self.wopfile_path = os.path.join(tempdir, 'Wopfile_{}.yml'.format(self.subcommand))
            (self.wopfile_path, wopfile_content) = self.create_wopfile(path=self.wopfile_path)
        wopmars_command_template = "wopmars -w {wopfile_path} -D sqlite:///{db} -p -v"
        wopmars_command = wopmars_command_template\
            .format(wopfile_path=self.wopfile_path, **self.parameters)
        if self.parameters['dryrun']:
            wopmars_command += " -n"
        if self.parameters['forceall']:
            wopmars_command += " -F"
        if not self.parameters['params'] is None:
            wopmars_command += " --params {params}".format(**self.parameters)
        if not self.parameters['targetrule'] is None:
            wopmars_command += " --targetrule {targetrule}".format(**self.parameters)
        ###################
        #
        # Wopmars command options of merge, otu and optimize
        #
        ###################
        # if self.subcommand == 'merge':
        #     wopmars_command += " --fastqinfo {fastqinfo} --fastqdir {fastqdir}"
        #     wopmars_command += " --fastainfo {fastainfo} --fastqdir {fastqdir}"
        # else:
        #     sys.exit(1)
        wopmars_command = wopmars_command.format(**self.parameters)
        return wopmars_command


    def merge(self):
        return "merge"