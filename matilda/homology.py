import collections
from matilda import helper_funcs
import math
import matildacpp


class PersistentHomologyComputer(object):
    """
    Class for computing persistent homology of a filtered chain complex or filtered simplicial complex.

    Attributes
    ----------
    simplices
        List containing all simplices in the complex, regardless of filtration value.
        Simplices are themselves lists of integers.

    dimension: int
        Dimension of largest simplex. Equal to length of largest simplex minus one.
        number_of_vertices: Number of vertices present in the complex, regardless of filtration value.

    simplices_indices
        Ordering of `simplices` consistent with filtration values of each simplex.

    appears_at
        Filtration values for all simplices in the complex.

    bars
        Dictionary indexed by homology degree, its values are dictionaries indexed by simplex indices, and their values
        are lists of two elements corresponding to start and end of a bar.
        Concisely:  ``bars[degree][simplex_index] = [start,end]``

    simplices_lookup
        Dictionary indexed by simplices with values their corresponding index in the filtration.



    References
    ----------
    .. [1] PAPER ON HOMOLOGY AND VECTOR ANNOTATIONS.
    .. [2] OTHER PAPER
    """

    def __init__(self):
        self.bars = collections.defaultdict(
            dict
        )  # BARS ARE PAIRS INDEXED BY SIMPLEX(INNER) AND DEGREE(OUTER)
        self.perscpp = None  # matildacpp.PersistentHomologyComputer()
        self.persistent_cycles = {}
        self.fsc = None

    def compute_persistent_homology(
        self, fsc, upper_degree=1, verbose=False, with_representatives=False, modulus=2
    ):
        self.fsc = fsc
        if modulus == 0:
            self.perscpp = matildacpp.PersistentHomologyComputerReal()
        else:
            self.perscpp = matildacpp.PersistentHomologyComputerMod()
        temp = matildacpp.FilteredSimplicialComplex()
        temp.appears_at = fsc.appears_at
        temp.dimension = fsc.dimension
        temp.simplices = fsc.simplices
        temp.simplices_indices = fsc.simplices_indices
        if with_representatives:
            self.perscpp.compute_persistent_homology(
                temp, upper_degree, verbose, modulus
            )
            self.persistent_cycles = self.perscpp.persistent_cycles
        else:
            self.perscpp.compute_persistent_homology_no_representatives(
                temp, upper_degree, verbose, modulus
            )
        self.bars = self.perscpp.bars
        self.proto_bars = self.perscpp.proto_bars

        self.boundary_matrix = self.perscpp.boundary_matrix
        self.reduced_boundary_matrix = self.perscpp.reduced_boundary_matrix

    def print_bars(self, dimension=None):
        """

        Parameters
        ----------
        dimension: int
            Rank of homology to print barcodes of (optional)

        Prints barcode

        Return
        ------
        void
        """

        def print_single_dimension_bars(dimension):
            if dimension in self.bars:
                print("Bars at dimension %d:" % dimension)
                if not self.bars[dimension]:
                    print("None")
                for y in self.bars[dimension]:
                    print(self.bars[dimension][y])
            else:
                print("Bars at dimension %d have not been computed" % dimension)

        if dimension is not None:
            if isinstance(dimension, int):
                dimension = [dimension]
        else:
            dimension = self.bars.keys()
        for x in dimension:
            print_single_dimension_bars(x)


def bottleneck_distance(barcode_1, barcode_2, max_value=None):
    """
    Computes the bottleneck_distance between two persistence diagrams (each passed as list-like collections of pairs)
    Parameters
    ----------
    barcode_1:
        List of pairs (a,b) comprising the bars in a barcode (or points in a persistence diagram)
    barcode_2:
        List of pairs (a,b) comprising the bars in a barcode (or points in a persistence diagram)
    max_value:
        Optional maximum value to deal with points at infinity. Defaults to the finite maximum value present in diagrams.
    Return
    ------
    Bottleneck distance between `barcode_1` and `barcode_2`.
    """

    import scipy.spatial.distance as ssd
    import scipy.optimize as so
    import numpy
    from typing import ValuesView

    if isinstance(barcode_1, dict):
        barcode_1 = list(barcode_1.values())
    if isinstance(barcode_2, dict):
        barcode_2 = list(barcode_2.values())
    if isinstance(barcode_1, ValuesView):
        barcode_1 = list(barcode_1)
    if isinstance(barcode_1, ValuesView):
        barcode_2 = list(barcode_2)
    barcode_1 = numpy.array(barcode_1)
    barcode_2 = numpy.array(barcode_2)
    if max_value is None:
        max_value = max(
            numpy.amax(barcode_1[numpy.isfinite(barcode_1)]),
            numpy.amax(barcode_2[numpy.isfinite(barcode_2)]),
        )
    if numpy.any(numpy.isinf(barcode_1)):
        print(
            "Warning: Infinite bars found. Clipping coordinate values to {}.".format(
                max_value
            )
        )
        barcode_1[barcode_1 == numpy.inf] = max_value
    if numpy.any(numpy.isinf(barcode_2)):
        print(
            "Warning: Infinite bars found. Clipping coordinate values to {}.".format(
                max_value
            )
        )
        barcode_2[barcode_2 == numpy.inf] = max_value
    distances = ssd.cdist(barcode_1, barcode_2)
    assignments = so.linear_sum_assignment(distances)
    result = distances[assignments].sum()
    return result
