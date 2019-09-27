# Create logger and set logger level from environment variable
import logging
import os
import sys


LOGGER_LEVEL = logging.NOTSET
LOGGER_LEVEL = logging.DEBUG
if "LOGGER_LEVEL" in os.environ:
    LOGGER_LEVEL = int(os.environ["LOGGER_LEVEL"])

logger = logging.getLogger("wopmetabarcoding")
logger.setLevel(LOGGER_LEVEL)
loggerHandlerStdout = logging.StreamHandler(sys.stdout)

WOPMETABARCODING_LOG = "wopmetabarcoding.log"
if "WOPMETABARCODING_LOG" in os.environ:
    WOPMETABARCODING_LOG = os.environ["WOPMETABARCODING_LOG"]
loggerHandlerFile = logging.FileHandler(WOPMETABARCODING_LOG)


formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
loggerHandlerStdout.setFormatter(formatter)
loggerHandlerFile.setFormatter(formatter)

logger.addHandler(loggerHandlerStdout)
logger.addHandler(loggerHandlerFile)

