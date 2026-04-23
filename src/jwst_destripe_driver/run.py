"""
Driver script for JWST destriping.

Usage::

    python run.py <config>

"""

import sys
from jwst_destripe_driver.destripe import destripe

destripe(sys.argv[1], verbose=True, testing=True)
