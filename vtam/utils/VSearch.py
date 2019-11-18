import subprocess

from vtam.utils.Logger import Logger

class VSearch():

    def __init__(self, parameters):
        """Creates a vsearch object based on the parameters

        :param parameters: Dictionary with key:value for store parameters or key:None for other parameters
        :return: void
        """
        self.parameters = parameters

    def run(self):
        """Run the vsearch

        :return: void
        """
        command = 'vsearch'
        for param in self.parameters:
            if not self.parameters[param] is None:
                command += ' {} {}'.format(param, self.parameters[param])
            else:
                command += ' {}'.format(param)
        Logger.instance().info(command)
        run_result = subprocess.run(command.split(), stdout=subprocess.PIPE)
        Logger.instance().info(run_result.stdout)

