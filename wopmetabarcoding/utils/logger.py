# Create logger and set logger level from environment variable
import logging
import os
import sys

LOGGER_LEVEL = logging.NOTSET
if "LOGGER_LEVEL" in os.environ:
    LOGGER_LEVEL = int(os.environ["LOGGER_LEVEL"])

logger = logging.getLogger("wopmetabarcoding")
logger.setLevel(LOGGER_LEVEL)
loggerHandlerStdout = logging.StreamHandler(sys.stdout)

loggerHandlerFile = logging.FileHandler('wopmetabarcoding.log')

formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
loggerHandlerStdout.setFormatter(formatter)
loggerHandlerFile.setFormatter(formatter)

logger.addHandler(loggerHandlerStdout)
logger.addHandler(loggerHandlerFile)

