import shlex
import subprocess

from vtam.utils.Logger import Logger


class VSearch(object):

    def __init__(self, parameters):
        """Creates a vsearch object based on the parameters

        :param parameters: Dictionary with key:value for store parameters or key:None for other parameters
        :return: void
        """
        self.parameters = parameters

    def create_command(self):
        """Create the vsearch command that will be run_name

        :return: void
        """

        command = 'vsearch'
        for param in self.parameters:
            if not self.parameters[param] is None:
                command += ' --{} {}'.format(param, self.parameters[param])
            else:
                command += ' --{}'.format(param)
        Logger.instance().debug(command)

        return command

    def run(self):
        """Run the vsearch

        :return: void
        """
        command = self.create_command()
        run_result = subprocess.run(shlex.split(command), capture_output=True)
        Logger.instance().info(run_result.stdout.decode())
        Logger.instance().info(run_result.stderr.decode())
