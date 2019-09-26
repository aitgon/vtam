# Create logger and set logger level from environment variable
import logging
import os
import sys

from wopmars.utils.ColorPrint import ColorPrint

from wopmetabarcoding.utils.OptionManager import OptionManager

LOGGER_LEVEL = logging.NOTSET
LOGGER_LEVEL = logging.DEBUG
if "LOGGER_LEVEL" in os.environ:
    LOGGER_LEVEL = int(os.environ["LOGGER_LEVEL"])

logger = logging.getLogger("VTAM")
logger.setLevel(LOGGER_LEVEL)

#####################
#
# Logger file
#
#####################
if "--log" in OptionManager.instance():
    s_path_log_file = OptionManager.instance()["--log"]
    loggerHandlerFile = logging.FileHandler(s_path_log_file)
    file_formatter = logging.Formatter('%(levelname)s :: %(asctime)s :: %(name)s :: %(message)s')
    loggerHandlerFile.setFormatter(file_formatter)
    logger.addHandler(loggerHandlerFile)

#####################
#
# Logger stdout
#
#####################
loggerHandlerStdout = logging.StreamHandler(sys.stdout)
stdout_formatter = logging.Formatter(ColorPrint.blue('%(levelname)s :: %(asctime)s :: %(name)s :: %(message)s'))
loggerHandlerStdout.setFormatter(stdout_formatter)
logger.addHandler(loggerHandlerStdout)

