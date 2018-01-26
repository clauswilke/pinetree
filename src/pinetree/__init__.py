from .pinetree import *

import os
import platform
import subprocess
import sys

PINETREE_DATA = os.path.join(os.path.dirname(__file__))

# Support running tests from the source tree
if not os.path.exists(PINETREE_DATA):
    _pinetree_data = os.path.abspath(os.path.join(
        os.path.dirname(__file__), '../_skbuild/cmake-install/pinetree'))
    if os.path.exists(_pinetree_data):
        PINETREE_DATA = _pinetree_data

PINETREE_BIN_DIR = os.path.join(PINETREE_DATA, 'bin')

def _program(name, args):
    return subprocess.call([os.path.join(PINETREE_BIN_DIR, name)] + args)

def pinetree_test():
    raise SystemExit(_program('pinetree_test', sys.argv[1:]))
