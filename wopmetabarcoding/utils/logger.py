# Create logger and set logger level from environment variable
import logging
import os
import sys

logger = logging.getLogger("wopmetabarcoding")
logger.setLevel(logging.NOTSET)
if "LOGGER_LEVEL" in os.environ:
    logger.setLevel(int(os.environ["LOGGER_LEVEL"]))
loggerHandlerStdout = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
loggerHandlerStdout.setFormatter(formatter)
logger.addHandler(loggerHandlerStdout)
