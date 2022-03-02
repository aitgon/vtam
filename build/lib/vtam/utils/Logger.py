import logging
import os
import sys
from logging.handlers import RotatingFileHandler
from vtam.utils.Singleton import Singleton

from termcolor import colored


class LessThanFilter(logging.Filter):
    def __init__(self, exclusive_maximum, name=""):
        super(LessThanFilter, self).__init__(name)
        self.max_level = exclusive_maximum

    def filter(self, record):
        # non-zero return means we log this message
        return 1 if record.levelno < self.max_level else 0


class LoggerArguments(dict, Singleton):
    """
    The LoggerArguments is a single dictionary with the log_verbosity and log values

    """

    def __init__(self, *args, **kwargs):
        """
        :param: dict_args: dictionnary containing the arguments passed to the command line.
        :return: void
        """
        super().__init__(*args, **kwargs)


class Logger(Singleton):
    """
    This class defines the vtam logger.

    To use it in a wopmars wrapper, we need to pass the log_verbosity to the wrapper.
    Then we get add this CLIargumentDict.instance()['log_verbosity'] = int(self.option("log_verbosity"))
    """

    def __init__(self):

        self.__logger = logging.getLogger('vtam')
        self.formatter_str = '%(asctime)s :: %(levelname)s :: %(name)s :: %(message)s'
        formatter = logging.Formatter(self.formatter_str)
        self.__logger.setLevel(logging.DEBUG)  # set root's level

        #######################################################################
        #
        # Get verbosity and log file tsv_path
        #
        #######################################################################

        if 'log_verbosity' in LoggerArguments.instance():
            verbosity = int(LoggerArguments.instance()['log_verbosity'])
        elif 'VTAM_LOG_VERBOSITY' in os.environ:
            verbosity = int(os.environ['VTAM_LOG_VERBOSITY'])
        else:
            verbosity = 0

        if 'log' in LoggerArguments.instance():
            log_file_path = LoggerArguments.instance()['log']
        elif 'VTAM_LOG_FILE' in os.environ:
            log_file_path = os.environ['VTAM_LOG_FILE']
        else:
            log_file_path = None

        #######################################################################
        #
        # Stream stderr
        #
        #######################################################################

        self.stream_handler_stderr = logging.StreamHandler(stream=sys.stderr)
        self.stream_handler_stderr.setFormatter(formatter)
        self.stream_handler_stderr.setLevel(logging.WARNING)
        self.__logger.addHandler(self.stream_handler_stderr)

        #######################################################################
        #
        # Stream stdout
        #
        #######################################################################

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

        if log_file_path is not None and not log_file_path == 'None':

            # log file in append mode of size 1 Mo and 1 backup
            # handler equivalent to stream_handler in term of logging level but
            # write in .log file
            self.__file_handler_stdout = RotatingFileHandler(
                log_file_path, 'a', 1000000, 1)
            self.__file_handler_stdout.setFormatter(formatter)
            if verbosity <= 0:
                self.__file_handler_stdout.setLevel(logging.WARNING)
            if verbosity == 1:
                self.__file_handler_stdout.setLevel(logging.INFO)
            elif verbosity >= 2:
                self.__file_handler_stdout.setLevel(logging.DEBUG)
            # self.__file_handler_stdout.addFilter(LessThanFilter(logging.WARNING))
            self.__logger.addHandler(self.__file_handler_stdout)

            # err file in append mode of size 1 Mo and 1 backup
            # this handler will write everything in the .err file.
            self.__file_handler_stderr = RotatingFileHandler(
                "{}.err".format(log_file_path.rsplit(".", 1)[0]), 'a', 1000000, 1)
            self.__file_handler_stderr.setFormatter(formatter)
            self.__file_handler_stderr.setLevel(logging.WARNING)
            self.__logger.addHandler(self.__file_handler_stderr)

    def debug(self, msg):
        formatter_stream = logging.Formatter(colored(self.formatter_str, 'blue', attrs=['bold']))
        self.stream_handler_stderr.setFormatter(formatter_stream)
        self.stream_handler_stdout.setFormatter(formatter_stream)
        self.__logger.debug(msg)

    def info(self, msg):
        formatter_stream = logging.Formatter(colored(self.formatter_str, 'green', attrs=['dark']))
        self.stream_handler_stderr.setFormatter(formatter_stream)
        self.stream_handler_stdout.setFormatter(formatter_stream)
        self.__logger.info(msg)

    def warning(self, msg):
        formatter_stream = logging.Formatter(colored(self.formatter_str, 'yellow', attrs=['bold']))
        self.stream_handler_stderr.setFormatter(formatter_stream)
        self.stream_handler_stdout.setFormatter(formatter_stream)
        self.__logger.warning(msg)

    def error(self, msg):
        formatter_stream = logging.Formatter(colored(self.formatter_str, 'red', attrs=[]))
        self.stream_handler_stderr.setFormatter(formatter_stream)
        self.stream_handler_stdout.setFormatter(formatter_stream)
        self.__logger.error(msg)

    def critical(self, msg):
        formatter_stream = logging.Formatter(colored(self.formatter_str, 'red', attrs=[]))
        self.stream_handler_stderr.setFormatter(formatter_stream)
        self.stream_handler_stdout.setFormatter(formatter_stream)
        self.__logger.critical(msg)
