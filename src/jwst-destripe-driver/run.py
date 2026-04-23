"""
Driver script for JWST destriping.

Usage::

    python run.py <config>


"""

import sys
from roman_hlis_l2_driver.destripe_interface.destripe import destripe_all_layers


destripe_all_layers(sys.argv[1], verbose=True)
