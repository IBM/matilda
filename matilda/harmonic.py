import numpy as np
import scipy, scipy.linalg
from matilda import PersistentHomologyComputer
from tqdm.auto import tqdm
from time import time


def project_cycle(
    cycle,
    orthom,
    verbose=False,
):
    """project the vector 'cycle' on the orthogonal complement of the space spanned by 'boundary'

    Parameters
    ----------
    cycle : numpy.array

    orthom : numpy.array
        orthogonal matrix whose columns are form an othonormal basis for the boundary

    verbose : bool, optional
        prints the dimension of the inputs, for debugging purposes, by default False

    Returns
    -------
    np.array
        the projected vector
    """

    if verbose:
        print(
            "projecting cycle of shape {} onto the orthogonal of the boundary {}".format(
                cycle.shape, orthom.shape
            )
        )

    harmonic = cycle.copy()

    for u in orthom.T:
        harmonic -= np.dot(cycle, u) * u

    return harmonic


def cycle_array_from_dict(cycle, fidx):
    """utility function to convert a cycle from a dict representation to an array

    Parameters
    ----------
    cycle : dict
        dict of the type {simplex_id: coefficient} representing a cycle as a linear combination of simplices
    fidx : list or numpy.array
        list of simplices' ids forming a basis for the space of cycles

    Returns
    -------
    np.array
        the vector representation of the cycle in the basis given by fidx
    """
    return np.array([cycle.get(face_id, 0) for face_id in fidx])


class HarmonicRepresentativesComputer:
    """class to compute the harmonic representatives cycles of a given barcode

    @misc{basu2022harmonic,
          title={Harmonic Persistent Homology},
          author={Saugata Basu and Nathanael Cox},
          year={2022},
          eprint={2105.15170},
          archivePrefix={arXiv},
          primaryClass={math.AT}
        }
    """

    def __init__(self, K, homology=None):
        """
        Parameters
        ----------
        K : matilda.prototyping.FilteredSimplicialComplex
            a filtered simplicial complex
        homology : matilda.PersistentHomologyComputer, optional
            the precomputed mod0 persistent homology of the complex, the parameter 'with_representatives' has to be set to 'True'
            by default None, if None the homology will be computed internally
        """

        self.K = K

        if homology is None:
            print("computing mod0 homology")
            self.homology = PersistentHomologyComputer()
            self.homology.compute_persistent_homology(
                K, with_representatives=True, modulus=0
            )
        else:
            self.homology = homology

        self.harmonic_cycles = dict()

    def compute_harmonic_cycles(
        self,
        dim,
        selected_cycles=None,
        reduced=True,
        precompute=False,
        verbose=0,
        **kwargs
    ):
        """compute the harmonic representative cycle for each cycle of dimension 'dim' and stores them in the nested dict 'self.harmonic_cycles'

        Parameters
        ----------
        dim : int
            dimension
        selected_cycles : list of ints
            ids of the cycles to project
        verbose : int, optional
            0, 1, 2, 3 for increasing levels of infos.
            by default 0, no info
        """
        self.harmonic_cycles[dim] = dict()

        if reduced:
            boundary_dict = self.homology.reduced_boundary_matrix
            if verbose >= 1:
                print("using the column reduced boundary matrix")
        else:
            boundary_dict = self.homology.boundary_matrix
            if verbose >= 1:
                print("using the full boundary matrix")

        if selected_cycles is None:
            cycles_ids = sorted(self.homology.persistent_cycles[dim].keys())
            if verbose >= 1:
                print(
                    "no cycles selected, projecting all {} {}-dimensional cycles".format(
                        len(cycles_ids), dim
                    )
                )
        else:
            # we need to add all other cycles that are still alive for each cycle in selected_cycles
            to_add = selected_cycles.copy()
            for z in selected_cycles:
                still_alive_for_z = [
                    id
                    for id in self.homology.persistent_cycles[dim]
                    if (
                        (self.homology.bars[dim][id][0] < self.homology.bars[dim][z][0])
                        and (
                            self.homology.bars[dim][id][1]
                            > self.homology.bars[dim][z][0]
                        )
                    )
                ]

                to_add += still_alive_for_z

            cycles_ids = sorted(list(set(to_add)))

            if verbose >= 1:
                print(
                    "{} cycles selected, we need to project {} {}-dimensional cycles".format(
                        len(selected_cycles), len(cycles_ids), dim
                    )
                )

        if precompute:
            if verbose >= 1:
                print("precomputing boundary matrix")
                start = time()
            # compute boundary matrix
            M, id_rows, id_cols = self.K.get_boundary_matrix_at_filtration(
                simplex_id=cycles_ids[-1],  # the latest id
                dim=dim + 1,
                boundary_dict=boundary_dict,
                mode="economic",
                verbose=verbose >= 2,
            )

            Q, _ = scipy.linalg.qr(M, mode="economic")

            if verbose >= 1:
                print(
                    "B shape {}. Q shape {}. {:.3f}s".format(
                        M.shape, Q.shape, time() - start
                    )
                )

        pbar = tqdm(cycles_ids, disable=not (verbose >= 1))
        first_cycle = True
        for this_id in pbar:
            if not precompute:
                pbar.set_description("cycle {}, computing boundary...".format(this_id))
                this_cycle = self.homology.persistent_cycles[dim][this_id]

                if verbose >= 2:
                    print("\ncurrent cycle {}".format(this_id))
                    start = time()

                if first_cycle:
                    # compute boundary matrix
                    M, fidx, cidx = self.K.get_boundary_matrix_at_filtration(
                        simplex_id=this_id,
                        dim=dim + 1,
                        boundary_dict=boundary_dict,
                        mode="economic",
                        verbose=verbose >= 2,
                    )

                    # compute orthonormal basis for the boundary
                    boundary = OrthonormalBasis(M)

                    first_cycle = False

                else:
                    # update the basis of the boundary
                    newB, fidx, new_cidx = self.K.get_boundary_matrix_at_filtration(
                        simplex_id=this_id,
                        dim=dim + 1,
                        boundary_dict=boundary_dict,
                        mode="economic",
                        cofaces_to_skip=cidx,
                        verbose=verbose >= 2,
                    )

                    cidx = np.concatenate((cidx, new_cidx))
                    # add rows of 0s to the basis
                    boundary.basis = np.vstack(
                        (
                            boundary.basis,
                            np.zeros(
                                shape=(
                                    len(fidx) - len(boundary.basis),
                                    boundary.basis.shape[1],
                                ),
                                dtype=boundary.basis.dtype,
                            ),
                        )
                    )
                    boundary.extend_basis(newB, verbose=(verbose >= 3), **kwargs)

                pbar.set_description(
                    "cycle {}, boundary shape {}".format(this_id, boundary.basis.shape)
                )

                if verbose >= 2:
                    print(
                        "the orthonormalized boundary matrix has shape {}".format(
                            boundary.basis.shape
                        )
                    )

                if verbose >= 2:
                    print("{:.3f}s".format(time() - start))
                    print("projecting on the boundary...")
                    start = time()

                # get the this cycle array and project it to the boundary
                this_cycle_arr = cycle_array_from_dict(this_cycle, fidx).astype(float)

                if (len(cidx) > 0) and (len(fidx) > 0):
                    this_harmonic = project_cycle(
                        this_cycle_arr,
                        boundary.basis,
                        verbose=(verbose >= 3),
                    )
                else:
                    this_harmonic = this_cycle_arr

                if verbose >= 2:
                    print("{:.3f}s".format(time() - start))

            else:
                # the full basis for the boundary matrix is precomputed
                this_cycle = self.homology.persistent_cycles[dim][this_id]

                fidx = [r for r in id_rows if r <= this_id]
                cidx = [c for c in id_cols if c <= this_id]

                pbar.set_description(
                    "cycle {}, boundary shape ({}, {})".format(
                        this_id, len(fidx), len(cidx)
                    )
                )

                # get the this cycle array and project it to the boundary
                this_cycle_arr = cycle_array_from_dict(this_cycle, fidx).astype(float)

                if (len(cidx) > 0) and (len(fidx) > 0):
                    this_harmonic = project_cycle(
                        this_cycle_arr,
                        Q[: len(fidx), : len(cidx)],
                        verbose=(verbose >= 3),
                    )
                else:
                    this_harmonic = this_cycle_arr

            # now we need to project this to the orthogonal of other harmonic representatives that might be still alive
            # loop through all of the already computer harmonic cycles
            # and see if some of them are still alive

            if verbose >= 2:
                print("projecting on the other harmonic cycles...")
                start = time()

            still_alive = [
                id
                for id in self.harmonic_cycles[dim]
                if (
                    (
                        self.homology.bars[dim][id][0]
                        < self.homology.bars[dim][this_id][0]
                    )
                    and (
                        self.homology.bars[dim][id][1]
                        > self.homology.bars[dim][this_id][0]
                    )
                )
            ]

            if len(still_alive) > 0:
                other_harmonics = np.array(
                    [
                        cycle_array_from_dict(self.harmonic_cycles[dim][id], fidx)
                        for id in still_alive
                    ]
                ).T

                # normalize the columns of other_harmonics and use them as an orthonormal basis
                this_harmonic = project_cycle(
                    this_harmonic,
                    other_harmonics / np.linalg.norm(other_harmonics, axis=0),
                    verbose=(verbose >= 3),
                )

            if verbose >= 2:
                print("{:.3f}s".format(time() - start))

            # save the new harmonic to the dict
            self.harmonic_cycles[dim][this_id] = dict(zip(fidx, this_harmonic))


class OrthonormalBasis:
    """Class to compute and update an orthonormal basis given a set of vectors, as columns in a numpy array."""

    def __init__(self, B=None, check=False, atol=1e-10, **kwargs):
        """_summary_

        Parameters
        ----------
        B : numpy.array, optional
            a 2d array whose columns are the initial set of vectors, to be orthonormalized
            by default None
        check : bool, optional
            to check wheter B is already orthonormal, can take a bit of time
            by default False
        atol : float, optional
            the atol parameter for numpy.isclose() method, by default 1e-10
        """
        self.atol = atol

        if B is not None:
            if check:
                if self.is_orthonormal(B):
                    self.basis = B
                else:
                    self.basis = self.orthonormalize(B, **kwargs)
            else:
                self.basis = self.orthonormalize(B, **kwargs)

            self.clean_zeros(self.basis)

    def clean_zeros(self, A):
        """sets to 0 elements that are close to 0

        Parameters
        ----------
        A : numpy.array
            input
        """
        A[np.isclose(A, 0, atol=self.atol)] = 0

    def is_orthonormal(self, A):
        """check whether the input matrix is orthonormal

        Parameters
        ----------
        A : numpy.array
            input 2d array

        Returns
        -------
        bool
            True if A.T @ A is (close to) the identity
        """
        return np.isclose(A.T @ A, np.eye(A.shape[1])).all()

    def orthonormalize(self, A, method="qr", verbose=False):
        """given a set of vectors, as columns in a numpy.array A, computes an orthonormal basis for the space spanned by those vectors

        Parameters
        ----------
        A : numpy.array
            2d array of shape (MxN) whose columns describe the set of vectors
        method : {"qr", "svd"}, optional
            wheter to use `scipy.linalg.qr` or `scipy.linalg.svd`, by default "qr"
        verbose : bool, optional
            print the matricies' shapes, by default False

        Returns
        -------
        numpy.array
            an orthonormal matrix of shape (MxK) whith K=range(A)

        Raises
        ------
        ValueError
            if method not it {"qr", "svd"}
        """

        if len(A.shape) < 2:
            return np.atleast_2d(A / np.linalg.norm(A)).T

        if min(A.shape) == 0:
            return A

        if method == "qr":
            # QR factorization with pivoting
            Q, R, P = scipy.linalg.qr(
                A, mode="economic", pivoting=True, check_finite=False
            )

            if verbose:
                print("Q {} , R {}".format(Q.shape, R.shape))

            # the first k columns are the basis we are looking for
            slice = False
            for k in range(min(R.shape)):
                if np.isclose(R[k, k], 0):
                    slice = True
                    break

            # this is to distinguish the case in wich the last element of the diagonal
            # is 0 or not 0
            if slice:
                return Q[:, :k]
            else:
                return Q[:, : k + 1]

        elif method == "svd":
            return scipy.linalg.orth(A)

        else:
            raise ValueError("unkown method, use 'qr' or 'svd'")

    def extend_basis(self, toadd, extend_method="gs", verbose=False):
        """given a new set of vectors, add them to the current orthonormal basis `self.basis` if any of them is not in the current span

        Parameters
        ----------
        toadd : numpy.array
            2d array whose columns describe the set of vectors to add to the orthonormal basis
        check : bool, optional
            to check wheter `toadd` is already orthonormal, can take a bit of time
            by default False
        extend_method : {"gs", "qr"}, optional
            which method to use to exted the basis
            "gs" is Modified Gram-Smith
            "qr" computes the QR factorization of the matrix (self.basis, toadd) obtained by addind the new vectors to the right of the current basis
            by default "gs"
        verbose : bool, optional
            prints some information, by default False

        Raises
        ------
        ValueError
            _description_
        """
        if len(toadd.shape) < 2:
            toadd = np.atleast_2d(toadd / np.linalg.norm(toadd)).T

        if verbose:
            print(
                "adding {} vectors to the basis {}".format(
                    toadd.shape, self.basis.shape
                )
            )

        if extend_method == "gs":
            toadd = self.orthonormalize(toadd)

            if verbose:
                print("using gram-smith")
            for u in toadd.T:
                for b in self.basis.T:
                    u -= np.inner(b, u) * b

                if not np.allclose(u, 0):
                    self.basis = np.hstack(
                        (self.basis, np.atleast_2d(u / np.linalg.norm(u)).T)
                    )
            self.clean_zeros(self.basis)

        elif extend_method == "qr":
            if verbose:
                print("using QR")
            self.basis = self.orthonormalize(np.hstack((self.basis, toadd)))
        else:
            raise ValueError("extend_method not recognised, use one of {'gs', 'qr'}")
