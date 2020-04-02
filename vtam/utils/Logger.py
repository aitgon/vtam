import logging
import os
import sys
from logging.handlers import RotatingFileHandler

from vtam.utils.Singleton import Singleton
from vtam.utils.OptionManager import OptionManager

from termcolor import colored

class LessThanFilter(logging.Filter):
    def __init__(self, exclusive_maximum, name=""):
        super(LessThanFilter, self).__init__(name)
        self.max_level = exclusive_maximum

    def filter(self, record):
        #non-zero return means we log this message
        return 1 if record.levelno < self.max_level else 0

class Logger(Singleton):
    """
    This class defines the vtam logger.

    To use it in a wopmars wrapper, we need to pass the log_verbosity to the wrapper.
    Then we get add this OptionManager.instance()['log_verbosity'] = int(self.option("log_verbosity"))
    """

    def __init__(self):

        self.__logger = logging.getLogger('vtam')
        self.formatter_str = '%(asctime)s :: %(levelname)s :: %(name)s :: %(message)s'
        formatter = logging.Formatter(self.formatter_str)
        self.__logger.setLevel(logging.DEBUG)  # set root's level

        if 'log_verbosity' in OptionManager.instance():
            verbosity = int(OptionManager.instance()['log_verbosity'])
        elif not os.getenv('VTAM_LOG_VERBOSITY') is None:
            verbosity = int(os.getenv('VTAM_LOG_VERBOSITY'))
        else:
            try:
                # Get log_file from wopmars option manager
                from wopmars.utils.OptionManager import OptionManager as wopmars_option_manager
                verbosity = int(wopmars_option_manager.instance()['-v'])
            except KeyError:
                verbosity = 0

        ################################################################################################################
        #
        # Stream stderr
        #
        ################################################################################################################

        self.stream_handler_stderr = logging.StreamHandler(stream=sys.stderr)
        self.stream_handler_stderr.setFormatter(formatter)
        self.stream_handler_stderr.setLevel(logging.WARNING)
        self.__logger.addHandler(self.stream_handler_stderr)

        ################################################################################################################
        #
        # Stream stdout
        #
        ################################################################################################################

        self.stream_handler_stdout = logging.StreamHandler(stream=sys.stdout)
        self.stream_handler_stdout.setFormatter(formatter)
        if verbosity <= 0:
            self.stream_handler_stdout.setLevel(logging.WARNING)
        if verbosity == 1:
            self.stream_handler_stdout.setLevel(logging.INFO)
        elif verbosity >= 2:
            self.stream_handler_stdout.setLevel(logging.DEBUG)
        self.stream_handler_stdout.addFilter(LessThanFilter(logging.WARNING))
        self.__logger.addHandler(self.stream_handler_stdout)

        ################################################################################################################
        #
        # File stderr
        #
        ################################################################################################################

        log_file_path = None

        if 'log_file' in OptionManager.instance():
            log_file_path = str(OptionManager.instance()['log_file'])
        else:
            try:
                # Get log_file from wopmars option manager
                from wopmars.utils.OptionManager import OptionManager as wopmars_option_manager
                log_file_path = str(wopmars_option_manager.instance()['--log'])
            except KeyError:
                log_file_path = None

        if not log_file_path is None and not log_file_path == 'None':

            log_stdout_path = log_file_path.rsplit(".", 1)[0]

            # log file in append mode of size 1 Mo and 1 backup
            # handler equivalent to stream_handler in term of logging level but write in .log file
            self.__file_handler_stdout = RotatingFileHandler(log_stdout_path + ".log", 'a', 1000000, 1)
            self.__file_handler_stdout.setFormatter(formatter)
            if verbosity <= 0:
                self.__file_handler_stdout.setLevel(logging.WARNING)
            if verbosity == 1:
                self.__file_handler_stdout.setLevel(logging.INFO)
            elif verbosity >= 2:
                self.__file_handler_stdout.setLevel(logging.DEBUG)
            self.__file_handler_stdout.addFilter(LessThanFilter(logging.WARNING))
            self.__logger.addHandler(self.__file_handler_stdout)

            # err file in append mode of size 1 Mo and 1 backup
            # this handler will write everything in the .err file.
            self.__file_handler_stderr = RotatingFileHandler(log_stdout_path + ".err", 'a', 1000000, 1)
            self.__file_handler_stderr.setFormatter(formatter)
            self.__file_handler_stderr.setLevel(logging.WARNING)
            self.__logger.addHandler(self.__file_handler_stderr)

    def debug(self, msg):
        formatter_stream = logging.Formatter(colored(self.formatter_str, 'cyan', attrs=['bold']))
        self.stream_handler_stderr.setFormatter(formatter_stream)
        self.stream_handler_stdout.setFormatter(formatter_stream)
        self.__logger.info(msg)

    def info(self, msg):
        formatter_stream = logging.Formatter(colored(self.formatter_str, 'blue', attrs=['bold']))
        self.stream_handler_stderr.setFormatter(formatter_stream)
        self.stream_handler_stdout.setFormatter(formatter_stream)
        self.__logger.info(msg)

    def warning(self, msg):
        formatter_stream = logging.Formatter(colored(self.formatter_str, 'magenta', attrs=['bold']))
        self.stream_handler_stderr.setFormatter(formatter_stream)
        self.stream_handler_stdout.setFormatter(formatter_stream)
        self.__logger.info(msg)

    def error(self, msg):
        formatter_stream = logging.Formatter(colored(self.formatter_str, 'red', attrs=['bold']))
        self.stream_handler_stderr.setFormatter(formatter_stream)
        self.stream_handler_stdout.setFormatter(formatter_stream)
        self.__logger.info(msg)

    def critical(self, msg):
        formatter_stream = logging.Formatter(colored(self.formatter_str, 'red', attrs=['bold', 'reverse']))
        self.stream_handler_stderr.setFormatter(formatter_stream)
        self.stream_handler_stdout.setFormatter(formatter_stream)
        self.__logger.info(msg)

