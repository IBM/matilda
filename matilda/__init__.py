import numpy

MAX_DIM = 1000  # GLOBAL MAXIMUM OF DIMENSION AN FSC CAN HAVE FOR THE WHOLE PACKAGE
boundary_coefs = numpy.array([(-1) ** i for i in range(MAX_DIM)], dtype="int")

from .filtered_simplicial_complexes import *
from .helper_funcs import *
from .homology import *
from .examples import *
from .plot import *
from .summaries import *
from .mapper import *