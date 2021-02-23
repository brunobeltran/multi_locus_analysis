"""Module for fitting multi-particle stuff."""
from .stats import *
from .dataframes import *
from . import plotting
from . import finite_window


from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
