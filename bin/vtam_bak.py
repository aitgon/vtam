#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

__author__ = "Aitor Gonzalez"

# If running from within source directory,
# add '../wopmars' to sys.path.
_libdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../')
if os.path.isfile(os.path.join(_libdir, 'wopmetabarcoding', '__init__.py')):
    sys.path.insert(0, _libdir)

from wopmetabarcoding import VTAM

if __name__ == "__main__":
    VTAM(sys.argv)

