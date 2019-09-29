import logging
from logging.handlers import RotatingFileHandler

from wopmars.utils.ColorPrint import ColorPrint

from wopmetabarcoding.utils.Singleton import Singleton
from wopmetabarcoding.utils.OptionManager import OptionManager

class Logger(Singleton):
    """
    This class defines the vtam logger.

    To use it in a wopmars wrapper, we need to pass the log_verbosity to the wrapper.
    Then we get add this OptionManager.instance()['log_verbosity'] = int(self.option("log_verbosity"))
    """

    def __init__(self):
        self.__logger = logging.getLogger("VTAM")
        self.__formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(name)s :: %(message)s')
        #####################
        #
        # Logger stdout
        #
        #####################
        self.__stream_handler = logging.StreamHandler()
        self.__stream_handler.setFormatter(self.__formatter)
        self.__logger.addHandler(self.__stream_handler)
        #####################
        #
        # Logger file
        #
        #####################
        # log file in append mode of size 1 Mo and 1 backup
        # handler equivalent to stream_handler in term of logging level but write to .log file
        if not OptionManager.instance()["log_file"] is None:
            log_file_path = OptionManager.instance()["log_file"]
            self.__file_handler = RotatingFileHandler(log_file_path, 'a', 1000000, 1)
            self.__file_handler.setFormatter(self.__formatter)
            self.__logger.addHandler(self.__file_handler)

        #####################
        #
        # Set logger level
        #
        #####################
        if OptionManager.instance()['log_verbosity'] <= 0:
            self.__logger.setLevel(logging.WARNING)
        elif OptionManager.instance()['log_verbosity'] == 1:
            self.__logger.setLevel(logging.INFO)
        else:
            self.__logger.setLevel(logging.DEBUG)

    def debug(self, msg):
        self.__logger.debug(msg)

    def info(self, msg):
        self.__logger.info(msg)

    def warning(self, msg):
        self.__logger.warning(msg)

    def error(self, msg):
        self.__logger.error(msg)
