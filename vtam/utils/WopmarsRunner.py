import os
import pathlib
import tempfile

import jinja2
import yaml

from vtam.utils.PathManager import PathManager
from vtam.utils.Singleton import Singleton
from vtam.utils import constants


class WopmarsRunner(Singleton):

    def __init__(self, command, cli_args_dic):
        """

        :param command: takes one of two values: filter or optimize
        :param cli_args_dic: dictionnary (CLIargumentDict.instance()) with command
        """
        self.command = command
        self.cli_args_and_numerical_params = {}
        self.cli_args_and_numerical_params.update(cli_args_dic)

        ################################################################################################################
        #
        # Load default numerical cli_args_dic and overwrite with custom cli_args_dic
        #
        ################################################################################################################

        num_params_default_dic = constants.get_dic_params_default()
        self.cli_args_and_numerical_params.update(num_params_default_dic)

        if not (cli_args_dic['params'] is None):
            with open(cli_args_dic['params']) as fin:
                num_params_user_dic = yaml.load(fin, Loader=yaml.SafeLoader)
                self.cli_args_and_numerical_params.update(num_params_user_dic)
        #
        self.wopfile_path = None
        self.tempdir = PathManager.instance().get_tempdir()


    def create_wopfile(self, path=None):
        """

        :param wopfile_path: Path of output wopfile
        :return: tuple (wopfile_path, wopfile_content)
        """

        ################################################################################################################
        #
        # Get Wopfile template output
        # Create Wopfile output
        #
        ################################################################################################################

        wopfile_path = path
        if wopfile_path is None:
            wopfile_path = tempfile.NamedTemporaryFile().name
        self.wopfile_path = wopfile_path

        ################################################################################################################
        #
        # Create Wopfile content
        #
        ################################################################################################################

        template_dir = os.path.join(os.path.dirname(__file__), '../data')
        jinja2_env = jinja2.Environment(loader=jinja2.FileSystemLoader(template_dir))
        template = None

        ################################################################################################################
        #
        # Filter command
        #
        ################################################################################################################

        if self.command == 'filter':
            template = jinja2_env.get_template('wopfile_filter.yml')
            # Create wopfile
            wopfile_path = os.path.join(self.tempdir, 'wopfile_filter.yml')

        ################################################################################################################
        #
        # Optimize command
        #
        ################################################################################################################

        elif self.command == 'optimize':
            outdir = self.cli_args_and_numerical_params['outdir']
            self.cli_args_and_numerical_params['sortreads'] = os.path.join(outdir, "sortreads.tsv")
            #
            self.cli_args_and_numerical_params['optimize_lfn_biosample_replicate'] \
                = os.path.join(outdir, "optimize_lfn_biosample_replicate.tsv")
            self.cli_args_and_numerical_params['optimize_lfn_read_count_and_lfn_variant'] \
                = os.path.join(outdir, "optimize_lfn_read_count_and_lfn_variant.tsv")
            self.cli_args_and_numerical_params['optimize_lfn_variant_specific'] \
                = os.path.join(outdir, "optimize_lfn_variant_specific.tsv")
            self.cli_args_and_numerical_params['optimize_pcr_error'] \
                = os.path.join(outdir, "optimize_pcr_error.tsv")
            self.cli_args_and_numerical_params['optimize_lfn_read_count_and_lfn_variant_replicate'] \
                = os.path.join(outdir, "optimize_lfn_read_count_and_lfn_variant_replicate.tsv")
            self.cli_args_and_numerical_params['optimize_lfn_variant_replicate_specific'] \
                = os.path.join(outdir, "optimize_lfn_variant_replicate_specific.tsv")
            #
            template = jinja2_env.get_template('wopfile_optimize.yml')
            wopfile_path = os.path.join(self.tempdir, 'wopfile_optimize.yml')

        wopfile_content = template.render(self.cli_args_and_numerical_params)

        ################################################################################################################
        #
        # Write to wopfile
        #
        ################################################################################################################

        pathlib.Path(os.path.dirname(wopfile_path)).mkdir(parents=True, exist_ok=True)
        with open(wopfile_path, "w") as fout:
            fout.write(wopfile_content)
        return wopfile_path, wopfile_content


    def get_wopmars_command(self):
        """
        Construct the wopmars command depending on the user arguments

        :return: string with output output of wopfile
        """

        ################################################################################################################
        #
        # Base wopmars command
        #
        ################################################################################################################

        if self.wopfile_path is None:
            self.wopfile_path = os.path.join(self.tempdir, 'Wopfile_{}.yml'.format(self.command))
            self.wopfile_path, wopfile_content = self.create_wopfile(path=self.wopfile_path)

        # Base command
        wopmars_command_template = "wopmars -w {wopfile_path} -D sqlite:///{db} "
        wopmars_command = wopmars_command_template\
            .format(wopfile_path=self.wopfile_path, **self.cli_args_and_numerical_params)

        if 'dryrun' in self.cli_args_and_numerical_params:
            if self.cli_args_and_numerical_params['dryrun']:
                wopmars_command += " -n"

        if 'forceall' in self.cli_args_and_numerical_params:
            if self.cli_args_and_numerical_params['forceall']:
                wopmars_command += " -F"

        if 'log_verbosity' in self.cli_args_and_numerical_params:
            if self.cli_args_and_numerical_params['log_verbosity'] > 0: # -v then pass this verbosity to wopmars
                wopmars_command += " -v"
                if self.cli_args_and_numerical_params['log_verbosity'] > 1: # -vv or higher, then do no pass it through environmental variables
                    os.environ['VTAM_LOG_VERBOSITY'] = str(self.cli_args_and_numerical_params['log_verbosity'])

        # if not 'log_file' in self.cli_args_and_numerical_params or self.cli_args_and_numerical_params['log_file'] is None:
        #     self.cli_args_and_numerical_params['log_file'] = '.vtam.log'
        if not self.cli_args_and_numerical_params['log_file'] is None:
            wopmars_command += " --log " + self.cli_args_and_numerical_params['log_file']

        if 'since' in self.cli_args_and_numerical_params:
            if not self.cli_args_and_numerical_params['since'] is None:
                wopmars_command += " --since {since}".format(**self.cli_args_and_numerical_params)

        if 'until' in self.cli_args_and_numerical_params:
            if not self.cli_args_and_numerical_params['until'] is None:
                wopmars_command += " --until {until}".format(**self.cli_args_and_numerical_params)

        wopmars_command = wopmars_command.format(**self.cli_args_and_numerical_params)
        return wopmars_command
