import logging

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
        self.__formatter = logging.Formatter(ColorPrint.blue('%(levelname)s :: %(asctime)s :: %(name)s :: %(message)s'))
        #
        self.__handler = logging.StreamHandler()
        self.__handler.setFormatter(self.__formatter)
        self.__logger.addHandler(self.__handler)
        self.seen = False

        #####################
        #
        # Set logger level
        #
        #####################
        if OptionManager.instance()['log_verbosity'] == 0:
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
