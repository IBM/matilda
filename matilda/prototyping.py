from .sparse_linear_algebra import SparseMatrixR, SparseVector, ModIntType
from copy import copy
import itertools
import numpy as np
from tqdm.auto import tqdm

"""
Returns a list with the simplex facets in the canonical order (to be used in computing in boundary).
"""


def simplex_facets(simplex):
    result = []
    for i, _ in enumerate(simplex):
        facet = copy(simplex)
        del facet[i]
        result.append(facet)
    return result


"""
Highly inneficient version of a filtered simplicial complex class. Used for prototyping only.
"""


class FilteredSimplicialComplex(object):
    """
    Class for modelling a filtered (finite) simplicial complex.

    Attributes
    ----------
    simplices
        List containing all simplices in the complex, regardless of filtration value.
        Simplices are themselves numpy arrays of integers.

    dimension: int
        Dimension of largest simplex. Equal to length of largest simplex minus one.
        number_of_vertices: Number of vertices present in the complex, regardless of filtration value.

    simplices_indices
        List of ints encoding a strict ordering of `simplices` according to filtration values of each simplex. It
        requires simplices to have larger indices than their faces. This means that to iterate over simplices
        in filtration order, we do so by iterating over simplices_indices, so simplices[simplices_indices[i]]
        is the ith simplex added to the filtration.

    appears_at
        List of all filtration values for all simplices in the complex. Also indexed by simplices_indices.

    """

    def __init__(
        self, dimension=0, simplices=None, simplices_indices=None, appears_at=None
    ):
        """
        Default initializer. Sets all main attributes of simplex to default values.

        Returns
        -------
        void
        """
        self.dimension = dimension
        if simplices_indices is None:
            self.simplices_indices = []
        else:
            self.simplices_indices = simplices_indices
        if appears_at is None:
            self.appears_at = []
        else:
            self.appears_at = appears_at
        if simplices is None:
            self.simplices = []
        else:
            self.simplices = simplices

    def get_boundary_matrix_at_filtration(
        self,
        boundary_dict,
        value=None,
        simplex_id=None,
        dim=1,
        mode="economic",
        cofaces_to_skip=[],
        verbose=False,
    ):
        """_summary_

        Parameters
        ----------
        value : _type_, optional
            _description_, by default None
        simplex_id : _type_, optional
            _description_, by default None
        dim : int, optional
            dimension of the boundary map (C_{dim} to C_{dim-1}), by default 1

        mode : {'full', 'economic'}, optional
            _description_, by default 'economic'

        verbose : bool, optional
            _description_, by default False

        Returns
        -------
        _type_
            _description_

        Raises
        ------
        ValueError
            _description_
        """
        if value:
            indices_list = self.subcomplex_at_filtration(value)
        elif simplex_id:
            indices_list = self.subcomplex_at_index(simplex_id)
        else:
            raise ValueError("specify either value or simplex_id")

        if mode == "economic":
            faces_idx = np.array(
                [
                    i
                    for i in indices_list
                    if len(self.simplices[self.simplices_indices[i]]) == dim
                ]
            )
        elif mode == "full":
            # returns all (d-1)-dimensional simplices ids
            faces_idx = np.array(
                [
                    i
                    for i in self.simplices_indices
                    if len(self.simplices[self.simplices_indices[i]]) == dim
                ]
            )
        else:
            raise ValueError("unknow mode, use one of {'full', 'economic'}")

        cofaces_idx = np.array(
            [
                i
                for i in indices_list
                if (len(self.simplices[self.simplices_indices[i]]) == dim + 1)
                and (i not in cofaces_to_skip)
                and (i in boundary_dict)
            ]
        )

        nrows = len(faces_idx)
        ncols = len(cofaces_idx)
        if verbose:
            print("creating {}x{} boundary matrix".format(nrows, ncols))

        boundary_matrix = np.zeros((nrows, ncols), dtype=int)

        for j, cid in enumerate(cofaces_idx):
            for fid in boundary_dict[cid]:
                i = np.where(faces_idx == fid)[0]
                boundary_matrix[i, j] = boundary_dict[cid][fid]

        return boundary_matrix, faces_idx, cofaces_idx

    def subcomplex_at_filtration(self, t):
        indices_list = []

        for i, id in enumerate(self.simplices_indices):
            if self.appears_at[id] <= t:
                indices_list.append(i)

        return indices_list

    def subcomplex_at_index(self, id):
        import warnings

        if id < 0:
            warnings.warn("{} is less than 0 - returning all simplices")
            return [i for i in self.simplices_indices]
        else:
            return [i for i in self.simplices_indices if i <= id]


"""
Class modelling a formal simplicial chain, distinct from a SparseVector.
"""


class SimplicialChain:
    def __init__(self, coefficients=None, simplices=None):
        if coefficients is None:
            coefficients = []
        if simplices is None:
            simplices = []
        self.coefficients = coefficients
        self.simplices = simplices

    @classmethod
    def from_pairs(cls, pairs):
        coefficients = []
        simplices = []
        for x, y in pairs:
            coefficients.append(x)
            simplices.append(y)
        return cls(coefficients, simplices)


class FilteredSimplicialChainComplex:
    def __init__(self, fsc, coeff_type=ModIntType(2)) -> None:
        self.fsc = fsc

    """
    Given a chain, it returns its algebraic boundary. As a formal operator it does not need to be associated
    with a Filtered Simplicial Complex, but for a matrix representation we need to fix an ordered basis.
    """

    def boundary(simplex, lookup, coeff_type=ModIntType(2)):
        coeffs = []
        terms = []
        accum = -1
        for i, _ in enumerate(simplex):
            facet = copy(simplex)
            del facet[i]
            coeffs.append(accum * (-1))
            terms.append(facet)
            accum = accum * (-1)
        return SparseVector.from_list(coeff_type, list(zip(terms, lookup[facet])))


class SimplexLookup:
    def __init__(self, n_points, max_dim):
        self.simplex_to_index = {
            v: k for k, v in enumerate(itertools.combinations(range(n_points), max_dim))
        }
        self.index_to_simplex = {k: v for v, k in self.simplex_to_index}
