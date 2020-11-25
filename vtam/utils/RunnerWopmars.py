import jinja2
import os
import pathlib
import tempfile

from vtam.utils.FileParams import FileParams
from vtam.utils.PathManager import PathManager
from vtam.utils.Singleton import Singleton


class RunnerWopmars(Singleton):

    def __init__(self, command, cli_args_dic):
        """

        :param command: takes one of two values: filter or optimize
        :param cli_args_dic: dictionnary (CLIargumentDict.instance()) with command
        """

        self.command = command
        self.cli_args_and_numerical_params = {}
        self.cli_args_and_numerical_params.update(cli_args_dic)

        # Add user params_lfn_variant.yml parameters

        params_dic = FileParams(cli_args_dic['params']).get_params_dic()
        self.cli_args_and_numerical_params.update(params_dic)

        self.wopfile_path = None
        self.tempdir = PathManager.instance().get_tempdir()

    def create_wopfile(self, path=None):
        """

        :param wopfile_path: Path of output wopfile
        :return: tuple (wopfile_path, wopfile_content)
        """

        ############################################################################################
        #
        # Create wopfile_path
        #
        ############################################################################################

        wopfile_path = path
        if wopfile_path is None:
            wopfile_path = tempfile.NamedTemporaryFile().name
        self.wopfile_path = wopfile_path

        ############################################################################################
        #
        # Load jinja2 environment
        #
        ############################################################################################

        template_dir = os.path.join(os.path.dirname(__file__), '../data')
        jinja2_env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(template_dir))
        template = None

        ############################################################################################
        #
        # Set either lfn_variant or lfn_variant_replicate algorithms
        #
        ############################################################################################

        if not self.cli_args_and_numerical_params['lfn_variant_replicate']:
            self.cli_args_and_numerical_params['lfn_variant_replicate_cutoff'] = None

        ############################################################################################
        #
        # Filter command
        #
        ############################################################################################

        if self.command == 'filter':
            template = jinja2_env.get_template('wopfile_filter.yml')
            # Create wopfile
            wopfile_path = os.path.join(self.tempdir, 'wopfile_filter.yml')

            if self.cli_args_and_numerical_params['lfn_variant_replicate']:  # lfn_variant replicate algorithm
                self.cli_args_and_numerical_params['lfn_variant_cutoff'] = None
                self.cli_args_and_numerical_params['lfn_variant_specific_cutoff'] = None
                if self.cli_args_and_numerical_params['cutoff_specific'] is None:
                    self.cli_args_and_numerical_params['lfn_variant_replicate_specific_cutoff'] = None
                else:
                    self.cli_args_and_numerical_params['lfn_variant_replicate_specific_cutoff'] = self.cli_args_and_numerical_params['cutoff_specific']


            else:  # lfn_variant algorithm
                self.cli_args_and_numerical_params['lfn_variant_replicate_cutoff'] = None
                self.cli_args_and_numerical_params['lfn_variant_replicate_specific_cutoff'] = None
                if self.cli_args_and_numerical_params['cutoff_specific'] is None:
                    self.cli_args_and_numerical_params['lfn_variant_specific_cutoff'] = None
                else:
                    self.cli_args_and_numerical_params['lfn_variant_specific_cutoff'] = self.cli_args_and_numerical_params['cutoff_specific']

        ############################################################################################
        #
        # Optimize command
        #
        ############################################################################################

        elif self.command == 'optimize':
            outdir = self.cli_args_and_numerical_params['outdir']
            # self.cli_args_and_numerical_params['sortreads'] = os.path.join(
            #     outdir, "sortreads.tsv")
            #
            self.cli_args_and_numerical_params['optimize_lfn_sample_replicate'] \
                = os.path.join(outdir, "optimize_lfn_sample_replicate.tsv")
            self.cli_args_and_numerical_params['optimize_pcr_error'] \
                = os.path.join(outdir, "optimize_pcr_error.tsv")

            if not self.cli_args_and_numerical_params['lfn_variant_replicate']:  # lfn_variant algorithm

                self.cli_args_and_numerical_params['optimize_lfn_read_count_and_lfn_variant'] \
                    = os.path.join(outdir, "optimize_lfn_read_count_and_lfn_variant.tsv")
                self.cli_args_and_numerical_params['optimize_lfn_variant_specific'] \
                    = os.path.join(outdir, "optimize_lfn_variant_specific.tsv")

            else:  # lfn_variant_replicate algorithm

                self.cli_args_and_numerical_params['optimize_lfn_read_count_and_lfn_variant'] \
                    = os.path.join(outdir, "optimize_lfn_read_count_and_lfn_variant_replicate.tsv")
                self.cli_args_and_numerical_params['optimize_lfn_variant_specific'] \
                    = os.path.join(outdir, "optimize_lfn_variant_replicate_specific.tsv")

            template = jinja2_env.get_template('wopfile_optimize.yml')
            wopfile_path = os.path.join(self.tempdir, 'wopfile_optimize.yml')

        wopfile_content = template.render(self.cli_args_and_numerical_params)

        ############################################################################################
        #
        # Write to wopfile
        #
        ############################################################################################

        pathlib.Path(
            os.path.dirname(wopfile_path)).mkdir(
            parents=True, exist_ok=True)
        with open(wopfile_path, "w") as fout:
            fout.write(wopfile_content)

        return wopfile_path, wopfile_content

    def get_wopmars_command(self):
        """
        Construct the wopmars command depending on the user arguments

        :return: string with output output of wopfile
        """

        ############################################################################################
        #
        # Create wopfile
        #
        ############################################################################################

        if self.wopfile_path is None:
            self.wopfile_path = os.path.join(
                self.tempdir, 'Wopfile_{}.yml'.format(
                    self.command))
            self.wopfile_path, wopfile_content = self.create_wopfile(
                path=self.wopfile_path)

        ############################################################################################
        #
        # Create wopmars command
        #
        ############################################################################################

        # Base command
        wopmars_command_template = "wopmars -w {wopfile_path} -D sqlite:///{db} "
        wopmars_command = wopmars_command_template .format(
            wopfile_path=self.wopfile_path, **self.cli_args_and_numerical_params)

        if 'dryrun' in self.cli_args_and_numerical_params:
            if self.cli_args_and_numerical_params['dryrun']:
                wopmars_command += " -n"

        if 'forceall' in self.cli_args_and_numerical_params:
            if self.cli_args_and_numerical_params['forceall']:
                wopmars_command += " -F"

        if 'log_verbosity' in self.cli_args_and_numerical_params:
            # -v then pass this verbosity to wopmars
            if self.cli_args_and_numerical_params['log_verbosity'] > 0:
                wopmars_command += " -v"
                # -vv or higher, then do no pass it through environmental variables
                if self.cli_args_and_numerical_params['log_verbosity'] > 1:
                    os.environ['VTAM_LOG_VERBOSITY'] = str(
                        self.cli_args_and_numerical_params['log_verbosity'])

        if not self.cli_args_and_numerical_params['log'] is None:
            wopmars_command += " --log " + \
                self.cli_args_and_numerical_params['log']

        if 'since' in self.cli_args_and_numerical_params:
            if not self.cli_args_and_numerical_params['since'] is None:
                wopmars_command += " --since {since}".format(
                    **self.cli_args_and_numerical_params)

        if 'until' in self.cli_args_and_numerical_params:
            if not self.cli_args_and_numerical_params['until'] is None:
                wopmars_command += " --until {until}".format(
                    **self.cli_args_and_numerical_params)

        wopmars_command = wopmars_command.format(
            **self.cli_args_and_numerical_params)
        return wopmars_command
