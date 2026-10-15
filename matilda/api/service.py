"""Harmonic persistent homology as a single computation with JSON-ready output.

This module wraps the existing matilda pipeline (Vietoris-Rips construction, persistent
homology over the reals and :class:`matilda.harmonic.HarmonicRepresentativesComputer`)
without changing it. It validates inputs and converts the results to plain Python types,
naming simplices by their vertices instead of matilda's internal filtration positions.
"""

import math
import time

import numpy as np
import scipy.spatial.distance

import matilda
import matilda.harmonic
import matilda.prototyping

METRICS = (
    "euclidean",
    "sqeuclidean",
    "cityblock",
    "chebyshev",
    "cosine",
    "correlation",
    "braycurtis",
    "canberra",
)


def validate_inputs(
    points=None,
    distance_matrix=None,
    filtration=None,
    metric="euclidean",
    max_dimension=2,
    upper_bound=None,
    homology_dimension=1,
    vertex_labels=None,
    max_points=None,
):
    """check the inputs of :func:`compute_harmonic_persistence` and normalize them

    Exactly one of `points`, `distance_matrix` or `filtration` must be given.

    Parameters
    ----------
    points : array-like, optional
        (n, d) point cloud; a Vietoris-Rips complex is built from the pairwise distances
    distance_matrix : array-like, optional
        (n, n) symmetric, non-negative matrix with zero diagonal; a Vietoris-Rips complex
        is built from it
    filtration : list of (simplex, value) pairs, optional
        an explicit filtered simplicial complex. Every face of every simplex must be listed
        with a value not larger than the simplex's value, and vertices must be 0..n-1.
        Simplices are ordered by value, then dimension, then their position in the list
    metric : str, optional
        distance used for `points`, one of `METRICS`, by default "euclidean"
    max_dimension : int, optional
        maximum simplex dimension of the Vietoris-Rips complex, by default 2.
        Ignored for `filtration`
    upper_bound : float, optional
        maximum edge length of the Vietoris-Rips complex, required for `points` and
        `distance_matrix`
    homology_dimension : int, optional
        dimension of the bars whose harmonic representatives are computed, by default 1
    vertex_labels : list, optional
        one label per vertex, used to name vertices in the output
    max_points : int, optional
        reject inputs with more vertices than this, by default no limit

    Returns
    -------
    dict
        the same keyword arguments, with inputs converted to numpy arrays (or, for
        `filtration`, to a list of (sorted vertex tuple, value) pairs in filtration order)

    Raises
    ------
    ValueError
        if the inputs are inconsistent
    """
    given = [
        name
        for name, value in (
            ("points", points),
            ("distance_matrix", distance_matrix),
            ("filtration", filtration),
        )
        if value is not None
    ]
    if len(given) != 1:
        raise ValueError(
            "give exactly one of 'points', 'distance_matrix' or 'filtration', got {}".format(
                given or "none"
            )
        )

    if isinstance(homology_dimension, bool) or not isinstance(
        homology_dimension, (int, np.integer)
    ):
        raise ValueError("'homology_dimension' must be an integer")
    if homology_dimension < 0:
        raise ValueError("'homology_dimension' must be non-negative")

    if filtration is not None:
        if upper_bound is not None:
            raise ValueError("'upper_bound' only applies to Vietoris-Rips inputs")
        filtration, n_vertices, dimension = _validate_filtration(filtration)
        if homology_dimension > dimension:
            raise ValueError(
                "'homology_dimension' {} is larger than the complex dimension {}".format(
                    homology_dimension, dimension
                )
            )
    else:
        if isinstance(max_dimension, bool) or not isinstance(
            max_dimension, (int, np.integer)
        ):
            raise ValueError("'max_dimension' must be an integer")
        if max_dimension < 1:
            raise ValueError("'max_dimension' must be at least 1")
        if homology_dimension >= max_dimension:
            # in a Vietoris-Rips complex truncated at max_dimension every top-dimensional
            # cycle is an infinite bar, so homology there says nothing about the data
            raise ValueError(
                "'homology_dimension' must be smaller than 'max_dimension' ({})".format(
                    max_dimension
                )
            )
        if upper_bound is None:
            raise ValueError(
                "'upper_bound' is required for 'points' and 'distance_matrix'"
            )
        if not (math.isfinite(upper_bound) and upper_bound > 0):
            raise ValueError("'upper_bound' must be a positive finite number")

        if points is not None:
            points = np.asarray(points, dtype=float)
            if points.ndim != 2 or points.shape[0] == 0 or points.shape[1] == 0:
                raise ValueError("'points' must be a non-empty 2d array")
            if not np.isfinite(points).all():
                raise ValueError("'points' contains non-finite values")
            if metric not in METRICS:
                raise ValueError(
                    "unknown metric '{}', use one of {}".format(metric, METRICS)
                )
            n_vertices = points.shape[0]
        else:
            distance_matrix = np.asarray(distance_matrix, dtype=float)
            if (
                distance_matrix.ndim != 2
                or distance_matrix.shape[0] != distance_matrix.shape[1]
                or distance_matrix.shape[0] == 0
            ):
                raise ValueError("'distance_matrix' must be a non-empty square matrix")
            if not np.isfinite(distance_matrix).all():
                raise ValueError("'distance_matrix' contains non-finite values")
            if (distance_matrix < 0).any():
                raise ValueError("'distance_matrix' contains negative values")
            if not np.allclose(distance_matrix, distance_matrix.T):
                raise ValueError("'distance_matrix' is not symmetric")
            if not np.allclose(np.diag(distance_matrix), 0):
                raise ValueError("'distance_matrix' must have a zero diagonal")
            n_vertices = distance_matrix.shape[0]

    if max_points is not None and n_vertices > max_points:
        raise ValueError(
            "the input has {} vertices, more than the limit of {}".format(
                n_vertices, max_points
            )
        )

    if vertex_labels is not None:
        vertex_labels = list(vertex_labels)
        if len(vertex_labels) != n_vertices:
            raise ValueError(
                "got {} vertex labels for {} vertices".format(
                    len(vertex_labels), n_vertices
                )
            )

    return dict(
        points=points,
        distance_matrix=distance_matrix,
        filtration=filtration,
        metric=metric,
        max_dimension=max_dimension,
        upper_bound=upper_bound,
        homology_dimension=homology_dimension,
        vertex_labels=vertex_labels,
        max_points=max_points,
    )


def _validate_filtration(filtration):
    """check an explicit filtration and return it sorted, with its vertex count and dimension"""
    values = {}
    for entry in filtration:
        try:
            simplex, value = entry
        except (TypeError, ValueError):
            raise ValueError("'filtration' entries must be (simplex, value) pairs")
        simplex = tuple(sorted(int(v) for v in simplex))
        if len(simplex) == 0:
            raise ValueError("'filtration' contains an empty simplex")
        if len(set(simplex)) != len(simplex):
            raise ValueError("simplex {} repeats a vertex".format(list(simplex)))
        if simplex[0] < 0:
            raise ValueError("simplex {} has a negative vertex".format(list(simplex)))
        value = float(value)
        if not math.isfinite(value):
            raise ValueError("simplex {} has a non-finite value".format(list(simplex)))
        if simplex in values:
            raise ValueError("simplex {} is listed twice".format(list(simplex)))
        values[simplex] = value

    if not values:
        raise ValueError("'filtration' is empty")

    for simplex, value in values.items():
        if len(simplex) == 1:
            continue
        for i in range(len(simplex)):
            face = simplex[:i] + simplex[i + 1 :]
            if face not in values:
                raise ValueError(
                    "face {} of simplex {} is missing".format(list(face), list(simplex))
                )
            if values[face] > value:
                raise ValueError(
                    "face {} appears after its simplex {}".format(
                        list(face), list(simplex)
                    )
                )

    vertices = sorted(s[0] for s in values if len(s) == 1)
    if vertices != list(range(len(vertices))):
        raise ValueError("vertices must be numbered 0..n-1 without gaps")

    # sorting by value, then dimension, puts every face before its cofaces; the sort
    # is stable, so ties keep the order of the list
    ordered = sorted(values.items(), key=lambda item: (item[1], len(item[0])))
    dimension = max(len(s) for s in values) - 1
    return ordered, len(vertices), dimension


def compute_harmonic_persistence(
    points=None,
    distance_matrix=None,
    filtration=None,
    metric="euclidean",
    max_dimension=2,
    upper_bound=None,
    homology_dimension=1,
    min_persistence=0.0,
    precompute=True,
    tolerance=1e-12,
    vertex_labels=None,
    max_points=None,
    max_simplices=None,
):
    """compute the barcode and the harmonic representatives of a filtered complex

    The complex is either a Vietoris-Rips complex (from `points` or `distance_matrix`) or
    an explicit `filtration`; see :func:`validate_inputs` for these and the other input
    parameters. Persistent homology is computed over the reals and harmonic representatives
    with :class:`matilda.harmonic.HarmonicRepresentativesComputer`.

    Parameters
    ----------
    min_persistence : float, optional
        only return harmonic representatives of bars at least this long, by default 0.
        Infinite bars are always returned. Bars still alive at the birth of a selected bar
        are also projected internally, as matilda requires
    precompute : bool, optional
        passed to `compute_harmonic_cycles`: orthonormalize the whole boundary matrix
        once instead of extending a basis bar by bar, by default True
    tolerance : float, optional
        drop chain coefficients with absolute value at most this, by default 1e-12.
        Vertex weights use all coefficients
    max_simplices : int, optional
        fail if the complex has more simplices than this, by default no limit

    Returns
    -------
    dict
        JSON-serializable results with keys

        - "input": a summary of the input and parameters
        - "n_simplices": number of simplices of each dimension
        - "vertex_labels": the given labels, or None
        - "barcode": every bar up to the meaningful dimension, as dicts with an opaque
          "id" (unique within the result), "dimension", "birth", "death" and
          "persistence" (None for infinite bars), and the "birth_simplex" and
          "death_simplex" that create and kill it
        - "harmonic": harmonic representatives of the selected bars, each with a
          "chain" of {"simplex", "coefficient"} terms
        - "weights": for each selected bar, the weight of every vertex, i.e. the sum of
          the absolute coefficients of the simplices containing it, as in Gurnari et al.,
          "Probing omics data via harmonic persistent homology"
        - "timings": seconds spent in each step

    Raises
    ------
    ValueError
        if the inputs are invalid, the complex exceeds `max_simplices`, or the harmonic
        computation produces non-finite coefficients
    """
    inputs = validate_inputs(
        points=points,
        distance_matrix=distance_matrix,
        filtration=filtration,
        metric=metric,
        max_dimension=max_dimension,
        upper_bound=upper_bound,
        homology_dimension=homology_dimension,
        vertex_labels=vertex_labels,
        max_points=max_points,
    )
    if not (math.isfinite(min_persistence) and min_persistence >= 0):
        raise ValueError("'min_persistence' must be a non-negative finite number")
    if not (math.isfinite(tolerance) and tolerance >= 0):
        raise ValueError("'tolerance' must be a non-negative finite number")
    if homology_dimension == 0 and not precompute:
        # compute_harmonic_cycles treats the bar with id 0, always present in degree 0,
        # as a missing argument when precompute is False
        raise ValueError("'precompute' must be true for 'homology_dimension' 0")

    timings = {}
    start = time.perf_counter()
    if inputs["filtration"] is not None:
        kind = "filtration"
        ordered = inputs["filtration"]
        n_vertices = sum(1 for s, _ in ordered if len(s) == 1)
        dimension = max(len(s) for s, _ in ordered) - 1
        top_bar_dimension = dimension
        K = matilda.prototyping.FilteredSimplicialComplex(
            dimension=dimension,
            simplices=[np.array(s) for s, _ in ordered],
            simplices_indices=list(range(len(ordered))),
            appears_at=[v for _, v in ordered],
        )
    else:
        if inputs["points"] is not None:
            kind = "points"
            dm = scipy.spatial.distance.squareform(
                scipy.spatial.distance.pdist(inputs["points"], metric=metric)
            )
            if not np.isfinite(dm).all():
                raise ValueError(
                    "metric '{}' produced non-finite distances".format(metric)
                )
        else:
            kind = "distance_matrix"
            dm = inputs["distance_matrix"]
        n_vertices = dm.shape[0]
        dimension = max_dimension
        top_bar_dimension = max_dimension - 1
        K_vr = matilda.FilteredSimplicialComplex()
        K_vr.construct_vietoris_from_metric(
            np.ascontiguousarray(dm, dtype=np.float32), max_dimension, upper_bound
        )
        # the harmonic computation needs the prototyping class
        K = matilda.prototyping.FilteredSimplicialComplex(
            dimension=dimension,
            simplices=K_vr.simplices,
            simplices_indices=K_vr.simplices_indices,
            appears_at=K_vr.appears_at,
        )
    timings["complex"] = time.perf_counter() - start

    n_simplices = [0] * (dimension + 1)
    for s in K.simplices:
        n_simplices[len(s) - 1] += 1
    if max_simplices is not None and len(K.simplices) > max_simplices:
        raise ValueError(
            "the complex has {} simplices, more than the limit of {}; "
            "lower 'upper_bound' or 'max_dimension'".format(
                len(K.simplices), max_simplices
            )
        )

    start = time.perf_counter()
    homology = matilda.PersistentHomologyComputer()
    homology.compute_persistent_homology(K, with_representatives=True, modulus=0)
    timings["persistent_homology"] = time.perf_counter() - start

    def vertices_of(position):
        return [int(v) for v in K.simplices[position]]

    def describe_bar(dim, bar_id):
        birth, death = homology.bars[dim][bar_id]
        finite = math.isfinite(death)
        return {
            "birth": float(birth),
            "death": float(death) if finite else None,
            "persistence": float(death - birth) if finite else None,
        }

    barcode = []
    for dim in sorted(homology.bars):
        if dim > top_bar_dimension:
            continue
        for bar_id in sorted(homology.bars[dim]):
            death_position = homology.proto_bars[dim][bar_id][1]
            barcode.append(
                {
                    "id": int(bar_id),
                    "dimension": int(dim),
                    **describe_bar(dim, bar_id),
                    "birth_simplex": vertices_of(bar_id),
                    "death_simplex": (
                        None if death_position == -1 else vertices_of(death_position)
                    ),
                }
            )

    H = homology_dimension
    cycles = homology.persistent_cycles.get(H, {})
    selected = [
        bar_id
        for bar_id in sorted(cycles)
        if homology.bars[H][bar_id][1] - homology.bars[H][bar_id][0] >= min_persistence
    ]

    start = time.perf_counter()
    harmonic_cycles = {}
    if selected:
        computer = matilda.harmonic.HarmonicRepresentativesComputer(K, homology)
        computer.compute_harmonic_cycles(
            dim=H,
            selected_cycles=None if len(selected) == len(cycles) else selected,
            precompute=precompute,
        )
        harmonic_cycles = computer.harmonic_cycles[H]
    timings["harmonic"] = time.perf_counter() - start

    harmonic = []
    weights = []
    for bar_id in selected:
        chain = harmonic_cycles[bar_id]
        if not all(math.isfinite(c) for c in chain.values()):
            raise ValueError(
                "the harmonic representative of bar {} has non-finite coefficients".format(
                    bar_id
                )
            )
        bar = {"bar_id": int(bar_id), **describe_bar(H, bar_id)}
        harmonic.append(
            {
                **bar,
                "chain": [
                    {"simplex": vertices_of(position), "coefficient": float(c)}
                    for position, c in sorted(chain.items())
                    if abs(c) > tolerance
                ],
            }
        )
        vertex_weights = np.zeros(n_vertices)
        for position, c in chain.items():
            vertex_weights[K.simplices[position]] += abs(c)
        weights.append({**bar, "weights": vertex_weights.tolist()})

    return {
        "input": {
            "kind": kind,
            "n_vertices": n_vertices,
            "dimension": dimension,
            "metric": metric if kind == "points" else None,
            "upper_bound": upper_bound,
            "homology_dimension": H,
            "min_persistence": min_persistence,
            "precompute": precompute,
        },
        "n_simplices": n_simplices,
        "vertex_labels": inputs["vertex_labels"],
        "barcode": barcode,
        "harmonic": harmonic,
        "weights": {
            "vertices": inputs["vertex_labels"] or list(range(n_vertices)),
            "bars": weights,
        },
        "timings": timings,
    }
