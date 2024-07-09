"""
Location of data files
======================

Use as ::

    from lipidorder.data.files import *

"""

__all__ = [
    "MEMBRANE_GRO",
    "MEMBRANE_XTC"
]

import importlib.resources

data_directory = importlib.resources.files("lipidorder") / "data"

MEMBRANE_GRO = data_directory / "membrane.gro"
MEMBRANE_XTC = data_directory / "membrane.xtc"
