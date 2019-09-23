#!/usr/bin/env python

from wopmetabarcoding.utils.ArgParser import ArgParser
from wopmetabarcoding.utils.OptionManager import OptionManager

from wopmetabarcoding.utils.OptionManager import OptionManager


class VTAM(object):

    def __init__(self, sys_argv):
        self.sys_argv = sys_argv
        parser = ArgParser.get_arg_parser()
        args = parser.parse_args(sys_argv)
        ##################################
        #
        # Store arguments in OptionsManager
        #
        ##################################
        for k in vars(args):
            OptionManager.instance()[k] = vars(args)[k]


if __name__ == '__main__':
    VTAM()
