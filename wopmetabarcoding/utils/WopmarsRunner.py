import tempfile

import jinja2
import os

import sys

from wopmetabarcoding.utils.Singleton import Singleton
from wopmetabarcoding.utils.constants import parameters_numerical
# from wopmetabarcoding.utils.utilities import tempdir


class WopmarsRunner(Singleton):

    def __init__(self, command, parameters):
        """

        :param command: takes one of three values: merge, otu or optimize
        :param parameters: dictionnary (OptionManager.instance()) with command
        """
        self.command = command
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
        #####################
        #
        # Get Wopfile template path
        # Create Wopfile path
        #
        #####################
        wopfile_path = path
        if wopfile_path is None:
            wopfile_path = tempfile.NamedTemporaryFile().name
        self.wopfile_path = wopfile_path
        # wopfile_template_path\
        #     = os.path.join(os.path.dirname(__file__), '../../data/Wopfile_{}.yml'.format(self.command))
        #####################
        #
        # Create Wopfile content
        #
        #####################
        template_dir = os.path.join(os.path.dirname(__file__), '../../data')
        jinja2_env = jinja2.Environment(loader=jinja2.FileSystemLoader(template_dir))
        if self.command == 'merge':
            template = jinja2_env.get_template('Wopfile_merge.yml')
        elif self.command in ['otu', 'optimize']:
            # Add path to sortreads file
            self.parameters['sortreads'] = os.path.abspath(os.path.join(self.parameters['outdir'], "sortreads.tsv"))
            if self.command == 'otu':
                self.parameters['otutable'] = os.path.abspath(os.path.join(self.parameters['outdir'], "otutable.tsv"))
                if self.parameters['threshold_specific']: # threshold variant specific
                    template = jinja2_env.get_template('Wopfile_otu_thresholdspecific.yml')
                else:
                    template = jinja2_env.get_template('Wopfile_otu_thresholdgeneral.yml')
                # Create wopfile
                wopfile_path = os.path.join(self.parameters['outdir'], 'Wopfile_otu.yml')
            elif self.command == 'optimize':
                #
                self.parameters['optimize_lfn_biosample_replicate'] \
                    = os.path.abspath(os.path.join(self.parameters['outdir'], "optimize_lfn_biosample_replicate.tsv"))
                self.parameters['optimize_lfn_read_count_and_lfn_variant'] \
                    = os.path.abspath(os.path.join(self.parameters['outdir'], "optimize_lfn_read_count_and_lfn_variant.tsv"))
                self.parameters['optimize_lfn_variant_specific'] \
                    = os.path.abspath(os.path.join(self.parameters['outdir'], "optimize_lfn_variant_specific.tsv"))
                self.parameters['optimize_pcr_error'] \
                    = os.path.abspath(os.path.join(self.parameters['outdir'], "optimize_pcr_error.tsv"))
                self.parameters['optimize_lfn_read_count_and_lfn_variant_replicate'] \
                    = os.path.abspath(os.path.join(self.parameters['outdir'], "optimize_lfn_read_count_and_lfn_variant_replicate.tsv"))
                self.parameters['optimize_lfn_variant_replicate_specific'] \
                    = os.path.abspath(os.path.join(self.parameters['outdir'], "optimize_lfn_variant_replicate_specific.tsv"))
                #
                template = jinja2_env.get_template('Wopfile_optimize.yml')
                wopfile_path = os.path.join(self.parameters['outdir'], 'Wopfile_optimize.yml')
        wopfile_content = template.render(self.parameters)
        ################
        #
        # Write to wopfile
        #
        ################
        with open(wopfile_path, "w") as fout:
            fout.write(wopfile_content)
        return wopfile_path, wopfile_content


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
            self.wopfile_path = os.path.join(tempdir, 'Wopfile_{}.yml'.format(self.command))
            (self.wopfile_path, wopfile_content) = self.create_wopfile(path=self.wopfile_path)
        wopmars_command_template = "wopmars -w {wopfile_path} -D sqlite:///{db} -p -v"
        wopmars_command = wopmars_command_template\
            .format(wopfile_path=self.wopfile_path, **self.parameters)
        if self.parameters['dryrun']:
            wopmars_command += " -n"
        if self.parameters['forceall']:
            wopmars_command += " -F"
        if not self.parameters['log_file'] is None:
            wopmars_command += " --log " + self.parameters['log_file']
        # if not self.parameters['params'] is None:
        #     wopmars_command += " --params {params}".format(**self.parameters)
        if not self.parameters['targetrule'] is None:
            wopmars_command += " --targetrule {targetrule}".format(**self.parameters)
        wopmars_command = wopmars_command.format(**self.parameters)
        return wopmars_command
