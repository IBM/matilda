import collections
import itertools
import copy
from dataclasses import dataclass
import numpy
from matilda.helper_funcs import param_parser
from abc import abstractmethod


class Number(type):
    def __init__(self, modulus):
        self.modulus = modulus

    @abstractmethod
    def __eq__(self, other):
        pass

    @abstractmethod
    def __ne__(self, other):
        pass

    @abstractmethod
    def __neg__(self):
        pass

    @abstractmethod
    def __add__(self, other):
        pass

    @abstractmethod
    def __mul__(self, other):
        pass

    @abstractmethod
    def __truediv__(self, other):
        pass

    @abstractmethod
    # Required only for multiplicative inverses (defined as 1/x, i.e. `other` should always be a unit or 1)
    def __rtruediv__(self, other):
        pass


class ModIntType:
    def __init__(self, modulus):
        self.modulus = modulus

    def __call__(self, value):
        return ModInt(value=value, modulus=self.modulus)

    def __repr__(self):
        return f"Modulo class {self.modulus}"


class ModInt:
    """
    Base class for modeling coefficients in a finite field. Different from `ModIntType` in that instancing a `ModIntType` object is a `ModInt` with the
    specified modulus. So for example with the following

    ```
    my_coefficients_mod_5 = ModIntType(5) 
    my_coefficients_mod_3 = ModIntType(3) 
    x = my_coefficients_mod_5(3)
    y = my_coefficients_mod_5(7)
    z = my_coefficients_mod_3(3)
    ```

    we have that `x` is 3 modulo 5, `y` is 7 modulo 5, and `z` is 3 modulo 3 (i.e. zero modulo 3).
    """

    def __init__(self, value, modulus):
        if isinstance(value, int):
            self.value = value
            self.modulus = modulus
            self.value = (self.value % self.modulus +
                          self.modulus) % self.modulus
        elif isinstance(value, ModInt):
            if modulus != value.modulus:
                raise ValueError(
                    f"Mismatched modulo equivalence classes: {modulus} and {value.modulus}.")
            self.value = value.value
            self.modulus = modulus

    def __eq__(self, other):
        if other == 0:  # We allow directly comparing to zero for expediency
            return self.value % self.modulus
        return (self.value - other.value) % self.modulus == 0

    def __ne__(self, other):
        return not self == other

    def __neg__(self):
        return ModInt(-self.value,
                      self.modulus)

    def __add__(self, other):
        return ModInt(self.value + other.value,
                      self.modulus)

    def __mul__(self, other):
        return ModInt(self.value * other.value,
                      self.modulus)

    def __truediv__(self, other):
        if other == ModInt(0, self.modulus):
            raise ZeroDivisionError()

        old_r, r = other.modulus, other.value
        old_t, t = 0, 1
        while r != 0:
            q = old_r // r
            old_t, t = t, old_t - q*t
            old_r, r = r, old_r - q*r
        if old_t < 0:
            old_t += self.modulus
        return ModInt(self.value * old_t, self.modulus)

    def __rtruediv__(self, other):
        if other == 1:
            return ModInt(1, self.modulus) / self

    def __repr__(self):
        return f"{self.value} mod {self.modulus}"


class SparseVector(object):
    """
    Base class for modelling sparse vectors with coefficients in a finite field.
    """

    def __init__(self, coeff_type=ModIntType(2)):
        # Keys are indices, values are values at indices
        self.entries = {}
        self.shape = [0]
        self.coeff_type = coeff_type

    def __getitem__(self, i):
        if i in self.entries.keys():
            return self.entries[i]
        else:
            return self.coeff_type(0)

    def __setitem__(self, i, x):
        x = self.coeff_type(x)
        if x != self.coeff_type(0):
            self.entries[i] = x
            self.shape[0] = max(self.shape[0], i+1)
        else:
            keys = self.entries.keys()
            if i in keys:
                del self.entries[i]
            self.shape[0] = max(self.entries.keys(), default=-1) + 1

    def __iter__(self):
        yield from self.entries.items()

    def __mul__(self, c):  # Scalar product
        result = SparseVector(self.coeff_type)
        for i, v in self:
            result[i] = v*c
        return result

    def __neg__(self):
        result = SparseVector(self.coeff_type)
        for k, v in self:
            result[k] = -v
        return result

    def __add__(self, b):
        result = SparseVector(self.coeff_type)
        all_keys = list(
            set.union(set(self.entries.keys()), set(b.entries.keys())))
        for k in all_keys:
            result[k] = self[k] + b[k]
        return result

    def __eq__(self, b):
        if set(self.entries.keys()) != set(b.entries.keys()):
            return False
        for k, v in self:
            if b[k] != v:
                return False
        return True

    @classmethod
    def from_list(cls, coeff_type, l):
        result = cls(coeff_type)
        for i, x in enumerate(l):
            result[i] = x
        return result

    @classmethod
    def from_dict(cls, coeff_type, l):
        result = cls(coeff_type)
        for k, v in l.items():
            result[k] = v
        return result

    def __str__(self):
        if self:
            to_print = numpy.zeros(self.shape, dtype="int")
            for k, v in self:
                to_print[k] = v.value
        return str(to_print)

    def __repr__(self):
        return self.__str__()

    def is_empty(self):
        if len(self.entries) == 0:
            return True
        return False
    
    def copy(self):
        return copy.deepcopy(self)

    def __len__(self):
        return len(self.entries)

class SparseMatrixR(object):
    """
    Class for modeling sparse matrices with compressed rows.
    Stores exactly all non-zero rows of a given matrix. Only maintains shape attribute.
    """

    def __init__(self, coeff_type=ModIntType(2)):
        self.rows = {}
        self.shape = [0, 0]
        self.coeff_type = coeff_type

    def __setitem__(self, index, item):
        if isinstance(index, tuple):
            i, j = index
            if i in self.rows.keys():
                self.rows[i][j] = item  # Item is a coefficient
            else:
                self.rows[i] = SparseVector.from_dict(
                    self.coeff_type, {j: item})
            if self.rows[i].is_empty():
                del self.rows[i]
        else:
            self.rows[index] = item  # Item is a vector
            if self.rows[index].is_empty():
                del self.rows[index]
        self.update_shape()

    def __getitem__(self, index):
        if isinstance(index, tuple):
            i, j = index
            if i in self.rows.keys():
                return self.rows[i][j]
        else:
            if index in self.rows.keys():
                return self.rows[index]
            else:
                return SparseVector.from_dict(self.coeff_type,{})

    def transpose(self):
        result = SparseMatrixR(coeff_type=self.coeff_type)
        for k, v in self.rows.items():
            for kk, vv in v.items():
                try:
                    result.rows[kk][k] = vv
                except:
                    result.rows[kk] = {k: vv}
        self.update_shape()
        return result

    def update_shape(self):
        if len(self.rows) == 0:
            self.shape = [0, 0]
        else:
            self.shape[0] = max(self.rows.keys()) + 1
            self.shape[1] = max([k for v in self.rows.values()
                                for (k, _) in v], default=-1) + 1

    def __str__(self):
        if self:
            self.update_shape()
            shape = self.shape
            to_print = numpy.zeros(shape, dtype="int")
            for k, v in self.rows.items():
                for kk, vv in v:
                    to_print[k][kk] = vv.value
            return str(to_print)
        else:
            return str(self)

    def __repr__(self):
        return self.__str__()

    def copy(self):
        return copy.deepcopy(self)

# Pivots are nonzero entries that are first according
# to some ordering on a given row.
# For example,
#
# 0 7 2 0
# 0 1 0 0
# 0 3 8 4
#
# The first row's principal entry according to `max` is
# 2, the second row's is 1, and the third one's is 4
# In `LinearSparseMatrixR` we record, for a given row `r`,
# its pivot as pivot[r] = (column, value)
# if a row doesn't have a pivot, we set pivot[r] = (-1, <undefined>)


class LinearSparseMatrixR(SparseMatrixR):
    """
    Class with methods tailored towards systems of linear equations through row reduction methods
    """

    def __init__(self, coeff_type=ModIntType(2), pivot_func = min):
        super().__init__(coeff_type=coeff_type)
        self.pivots = {}
        self.pivots_lookup = {}         # Lookup dict has form `{column : row}`, where 
                                        # `row` is such that `self[row, column]` is a pivot
                                        # (if there is more than one row satisfying this, the
                                        # choice is arbitrary)
        self.pivot_func = pivot_func 

    def __setitem__(self, index, item):
        super().__setitem__(index, item)
        self.update_pivots()

    def set_column(self, j, x): # Assigns SparseVector `x` to column `j`
        keys = set(self.rows.keys())
        for k, v in x:
            self[k, j] = v
            if k in keys:
                keys.remove(k)
        for k in keys:
            self[k, j] = 0
        
    def addscale_from_to(self, from_index, to_index, coefficient):
        self[to_index] = self[to_index] + \
            self[from_index] * coefficient

    def scale_row(self, row_index, coefficient):
        self[row_index] = self[row_index]*coefficient

    def update_pivot(self, j): # Computes and sets the pivot of row `j``
        column = self.pivot_func(self[j].entries.keys(), default = -1)
        self.pivots[j] = (column, self[j, column])

    def update_pivots(self):
        self.pivots = {}
        for k in self.rows.keys():
            self.update_pivot(k)
        self.update_pivot_lookup()

    def update_pivot_lookup(self): 
        self.pivots_lookup = {v:k for k,v in self.pivots.items()}

def gauss_jordan(m):
    """
    Performs *row-wise* Gauss-Jordan elimination on *m* without swapping rows and following
    the ordering on row elements.
    E.g. With *max* (pivots are rightmost nonzero entries) and

            1 1 1
    m   =   0 1 1
            1 0 0
            0 0 1

    we have
            1 1 1       1 1 1       0 0 1       0 0 1
    m ->    1 0 0   ->  1 0 0   ->  1 0 0   ->  1 0 0   =   gauss_jordan(m)
            1 0 0       1 0 0       1 0 0       0 0 0
            0 0 1       1 1 0       1 1 0       0 1 0
    """
    result = m.copy()





def gauss_jordan_no_swappingR(m, pivot_func=min):
    """
    Performs in-place modified Gauss Jordan elimination on *m* without swapping rows. 
    Thus it reduces a matrix to its row reduced form modulo row swapping. 
    Does not have return value unless specified. Optionally applies elementary row operations to
    a diagonal matrix with as many rows as *m* (this can only be used with in-place modifications for now).
    """
    result = m
    rows_remaining = set(result.rows.keys())
    result.update_pivots()
    while True:
        pivot_pair = pivot_func(
            ((k, result.pivots[k]) for k in rows_remaining), key=lambda x: x[1], default=(-1, numpy.inf))
        # pivot_pair = min(((k,v) for (k,v) in result.pivots.items() if k not in rows_reduced), key = lambda x: x[1], default= (-1,numpy.inf))
        if pivot_pair[0] == -1:
            break
        rows_remaining.remove(pivot_pair[0])
        keys_off_pivot = [(k, pivot_pair[1]) for (k, v) in result.rows.items(
        ) if k != pivot_pair[0] and pivot_pair[1] in v.keys()]
        result.scale_row(
            pivot_pair[0], modulus_inv_table[result.rows[pivot_pair[0]][pivot_pair[1]]])
        for k in keys_off_pivot:
            coeff = modulus_neg_table[result.rows[k[0]][k[1]]]
            result.addscale_from_to(pivot_pair[0], k[0], coeff)
            # if the whole row was cancelled we might as well consider it already reduced
            if k[0] not in result.rows.keys():
                rows_remaining.remove(k[0])
            else:
                result.update_pivot(k[0])


def matrix_multiplication(a, b):
    result = SparseMatrixR()
    a.update_shape()
    c = b.transpose()
    c.update_shape()
    for k, v in a.rows.items():
        result.rows[k] = {}
        for kk, vv in c.rows.items():
            entry = modulus_table[sum(v[i]*vv[i]
                                      for i in v.keys() if i in vv.keys())]
            if entry:
                result.rows[k][kk] = entry
        if not result.rows[k]:
            del result.rows[k]
    return result


def multiply_matrix_vector(A, b):
    result = SparseVector()
    for k, v in A.rows.items():
        accum = 0
        for kk, vv in v.items():
            accum = accum + b[kk]*vv
        if accum != 0:
            result[k] = modulus_table[accum]
    return result


def solveR(A, b, pivot_func=min):  # solves Ax = -b, returns empty sparsevector if no solutions exist
    if not A.rows:
        return SparseVector.from_dict({})
    A.update_shape()
    last_col = A.shape[1]
    A.set_column(last_col, b)
    A.update_shape()
    # print([v for k,v in A.rows.items()])
    gauss_jordan_no_swappingR(A, pivot_func=pivot_func)
    # print(A.shape)

    for i in A.rows.keys():
        if A.get_entry(i, last_col) != 0 and not [k for k in A.rows[i].keys() if k != last_col]:
            return SparseVector.from_dict({})
    result = SparseVector()
    for i in A.rows.keys():
        if A.get_entry(i, last_col) != 0:
            pair = max([(k, v) for (k, v) in A.rows[i].items() if k !=
                       last_col], key=lambda x: x[0], default=[-1, -1])
            result[pair[0]] = pair[1]
    return result


def gauss_jordan_no_swappingR_with_aug(m, aug, pivot_func=min):
    """
    Performs in-place modified Gauss Jordan elimination on *m* without swapping rows. 
    Thus it reduces a matrix to its row reduced form modulo row swapping. 
    Does not have return value unless specified. Optionally applies elementary row operations to
    a diagonal matrix with as many rows as *m* (this can only be used with in-place modifications for now).
    """
    result = m
    rows_reduced = []
    result.update_pivots()
    rows_remaining = sorted(result.rows.keys())
    while True:
        pivot_pair = pivot_func(
            ((k, result.pivots[k]) for k in rows_remaining), key=lambda x: x[1], default=(-1, numpy.inf))
        # pivot_pair = min(((k,v) for (k,v) in result.pivots.items() if k not in rows_reduced),
        #                 key = lambda x: x[1], default= (-1,numpy.inf))
        if pivot_pair[0] == -1:
            break
        rows_remaining.remove(pivot_pair[0])
        keys_off_pivot = [(k, pivot_pair[1]) for (k, v) in result.rows.items(
        ) if k != pivot_pair[0] and pivot_pair[1] in v.keys()]
        aug.scale_row(
            pivot_pair[0], modulus_inv_table[result.rows[pivot_pair[0]][pivot_pair[1]]])
        result.scale_row(
            pivot_pair[0], modulus_inv_table[result.rows[pivot_pair[0]][pivot_pair[1]]])
        for k in keys_off_pivot:
            coeff = modulus_neg_table[result.rows[k[0]][k[1]]]
            aug.addscale_from_to(pivot_pair[0], k[0], coeff)
            result.addscale_from_to(pivot_pair[0], k[0], coeff)
            if k[0] not in result.rows.keys():
                rows_remaining.remove(k[0])
            else:
                result.update_pivot(k[0])


class StatefulSolve(object):
    def __init__(self):
        self.previous = SparseMatrixR()
        self.transformation = SparseMatrixR()
        self.current_iteration = 0

    def solve(self, A, b):  # solves Ax = -b, returns empty sparsevector if no solutions exist
        if self.current_iteration == 0:
            self.previous = A.copy()
            self.transformation = SparseMatrixR()
        self.current_iteration += 1
        self.previous.update_shape()
        for k in range(self.transformation.shape[0], max(self.previous.shape[0], b.shape[0], A.shape[0])):
            self.transformation.rows[k] = {k: 1}
        self.transformation.update_shape()
        last_col = self.previous.shape[1]
        b = multiply_matrix_vector(self.transformation, b)
        self.previous.set_column(last_col, b)
        gauss_jordan_no_swappingR_with_aug(self.previous, self.transformation)
        for i in self.previous.rows.keys():
            if self.previous.get_entry(i, last_col) != 0 and not [k for k in self.previous.rows[i].keys() if k != last_col]:
                return SparseVector.from_dict({})
        result = SparseVector()
        for i in self.previous.rows.keys():
            if self.previous.get_entry(i, last_col) != 0:
                pair = max([(k, v) for (k, v) in self.previous.rows[i].items(
                ) if k != last_col], key=lambda x: x[0])
                result[pair[0]] = pair[1]
        self.previous.set_column(last_col, SparseVector())
        return result


def hconcat(a, b):
    result = SparseMatrixR()
    a.update_shape()
    for k, v in a.rows.items():
        result.rows[k] = {}
        for kk, vv in v.items():
            result.rows[k][kk] = vv
    for k, v in b.rows.items():
        if k not in result.rows.keys():
            result.rows[k] = {}
        for kk, vv in v.items():
            result.rows[k][kk+a.shape[1]] = vv
    return result


def gauss_no_swappingR(m, pivot_func=min):
    """
    Performs in-place modified Gauss Jordan elimination on *m* without swapping rows. 
    Thus it reduces a matrix to its row reduced form modulo row swapping. 
    Does not have return value unless specified. Optionally applies elementary row operations to
    a diagonal matrix with as many rows as *m* (this can only be used with in-place modifications for now).
    """
    result = m
    rows_remaining = set(result.rows.keys())
    result.update_pivots()
    while True:
        pivot_pair = pivot_func(
            ((k, result.pivots[k]) for k in rows_remaining), key=lambda x: x[1], default=(-1, numpy.inf))
        # pivot_pair = min(((k,v) for (k,v) in result.pivots.items() if k not in rows_reduced), key = lambda x: x[1], default= (-1,numpy.inf))
        if pivot_pair[0] == -1:
            break
        rows_remaining.remove(pivot_pair[0])
        keys_off_pivot = [(k, pivot_pair[1]) for (k, v) in result.rows.items(
        ) if k > pivot_pair[0] and pivot_pair[1] in v.keys()]
        result.scale_row(
            pivot_pair[0], modulus_inv_table[result.rows[pivot_pair[0]][pivot_pair[1]]])
        for k in keys_off_pivot:
            coeff = modulus_neg_table[result.rows[k[0]][k[1]]]
            result.addscale_from_to(pivot_pair[0], k[0], coeff)
            # if the whole row was cancelled we might as well consider it already reduced
            if k[0] not in result.rows.keys():
                rows_remaining.remove(k[0])
            else:
                result.update_pivot(k[0])
