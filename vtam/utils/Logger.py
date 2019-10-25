import logging
from logging.handlers import RotatingFileHandler

from vtam.utils.Singleton import Singleton
from vtam.utils.OptionManager import OptionManager

from termcolor import colored

class Logger(Singleton):
    """
    This class defines the vtam logger.

    To use it in a wopmars wrapper, we need to pass the log_verbosity to the wrapper.
    Then we get add this OptionManager.instance()['log_verbosity'] = int(self.option("log_verbosity"))
    """

    def __init__(self):
        self.__logger = logging.getLogger("VTAM")
        self.__formatter_str = '%(asctime)s :: %(levelname)s :: %(name)s :: %(message)s'
        self.__formatter = logging.Formatter(self.__formatter_str)
        #####################
        #
        # Logger stdout
        #
        #####################
        self.__stream_handler = logging.StreamHandler()
        self.__logger.addHandler(self.__stream_handler)
        #####################
        #
        # Logger file
        #
        #####################
        # log file in append mode of size 1 Mo and 1 backup
        # handler equivalent to stream_handler in term of logging level but write to .log file
        try:
            assert not (OptionManager.instance()["log_file"] is None or OptionManager.instance()["log_file"] == 'None')
            log_file_path = OptionManager.instance()["log_file"]
            self.__file_handler = RotatingFileHandler(log_file_path, 'a', 1000000, 1)
            self.__file_handler.setFormatter(self.__formatter)
            self.__logger.addHandler(self.__file_handler)
        except AssertionError: # ignore
            pass
        except KeyError: # ignore
            pass

        #####################
        #
        # Set logger level
        #
        #####################
        try:
            if OptionManager.instance()['log_verbosity'] <= 0:
                self.__logger.setLevel(logging.WARNING)
            elif OptionManager.instance()['log_verbosity'] == 1:
                self.__logger.setLevel(logging.INFO)
            else:
                self.__logger.setLevel(logging.DEBUG)
        except KeyError: # ignore
            pass

    def debug(self, msg):
        formatter_stream_str = colored(self.__formatter_str, 'yellow', attrs=['bold'])
        formatter_stream = logging.Formatter(formatter_stream_str)
        self.__stream_handler.setFormatter(formatter_stream)
        self.__logger.debug(msg)

    def info(self, msg):
        formatter_stream_str = colored(self.__formatter_str, 'blue', attrs=['bold'])
        formatter_stream = logging.Formatter(formatter_stream_str)
        self.__stream_handler.setFormatter(formatter_stream)
        self.__logger.info(msg)

    def warning(self, msg):
        formatter_stream_str = colored(self.__formatter_str, 'magenta', attrs=['bold'])
        formatter_stream = logging.Formatter(formatter_stream_str)
        self.__stream_handler.setFormatter(formatter_stream)
        self.__logger.warning(msg)

    def error(self, msg):
        formatter_stream_str = colored(self.__formatter_str, 'red', attrs=['bold'])
        formatter_stream = logging.Formatter(formatter_stream_str)
        self.__stream_handler.setFormatter(formatter_stream)
        self.__logger.error(msg)
