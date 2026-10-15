import itertools

import numpy as np
import pytest
from scipy.spatial import distance_matrix

import matilda
import matilda.harmonic
import matilda.prototyping
from matilda.api import compute_harmonic_persistence, validate_inputs

# the complex of examples/harmonic/example_harmonic_1.ipynb
EXAMPLE_SIMPLICES = [
    [0],
    [1],
    [2],
    [0, 1],
    [0, 2],
    [1, 2],
    [3],
    [0, 1, 2],
    [0, 3],
    [2, 3],
]
EXAMPLE_VALUES = [0, 1, 2, 3, 3, 3, 4, 5, 6, 6]
EXAMPLE_FILTRATION = list(zip(EXAMPLE_SIMPLICES, EXAMPLE_VALUES))

# boundary of the tetrahedron, a triangulated 2-sphere
SPHERE_FILTRATION = [
    (list(s), float(len(s)))
    for d in (1, 2, 3)
    for s in itertools.combinations(range(4), d)
]


def two_circles(n=15, seed=0):
    rng = np.random.default_rng(seed)
    t = rng.uniform(0, 2 * np.pi, size=2 * n)
    points = np.c_[np.cos(t), np.sin(t)] + 0.05 * rng.normal(size=(2 * n, 2))
    points[n:, 0] += 3
    return points


def direct_harmonic(K, H, precompute):
    """harmonic cycles computed with matilda directly, keyed by bar id and vertex tuple"""
    homology = matilda.PersistentHomologyComputer()
    homology.compute_persistent_homology(K, with_representatives=True, modulus=0)
    computer = matilda.harmonic.HarmonicRepresentativesComputer(K, homology)
    computer.compute_harmonic_cycles(dim=H, precompute=precompute)
    return {
        bar_id: {tuple(int(v) for v in K.simplices[i]): c for i, c in chain.items()}
        for bar_id, chain in computer.harmonic_cycles[H].items()
    }


def chains_by_bar(result):
    return {
        r["bar_id"]: {tuple(t["simplex"]): t["coefficient"] for t in r["chain"]}
        for r in result["harmonic"]
    }


def assert_same_chains(direct, served):
    assert direct.keys() == served.keys()
    for bar_id, chain in direct.items():
        kept = {s: c for s, c in chain.items() if c != 0}
        assert kept.keys() == served[bar_id].keys()
        for s, c in kept.items():
            assert served[bar_id][s] == pytest.approx(c)


@pytest.mark.parametrize("precompute", [False, True])
def test_filtration_matches_matilda(precompute):
    K = matilda.prototyping.FilteredSimplicialComplex(
        dimension=2,
        simplices=[np.array(s) for s in EXAMPLE_SIMPLICES],
        simplices_indices=list(range(len(EXAMPLE_SIMPLICES))),
        appears_at=EXAMPLE_VALUES,
    )
    result = compute_harmonic_persistence(
        filtration=EXAMPLE_FILTRATION, precompute=precompute, tolerance=0
    )
    assert_same_chains(direct_harmonic(K, 1, precompute), chains_by_bar(result))


def test_vietoris_rips_matches_matilda():
    points = two_circles()
    K_vr = matilda.FilteredSimplicialComplex()
    K_vr.construct_vietoris_from_metric(
        distance_matrix(points, points).astype(np.float32), 2, 1.5
    )
    K = matilda.prototyping.FilteredSimplicialComplex(
        dimension=2,
        simplices=K_vr.simplices,
        simplices_indices=K_vr.simplices_indices,
        appears_at=K_vr.appears_at,
    )
    result = compute_harmonic_persistence(points=points, upper_bound=1.5, tolerance=0)
    assert result["n_simplices"][0] == len(points)
    assert len(result["harmonic"]) >= 2  # one long bar per circle
    assert_same_chains(direct_harmonic(K, 1, True), chains_by_bar(result))


def test_points_and_distance_matrix_agree():
    points = two_circles()
    from_points = compute_harmonic_persistence(points=points, upper_bound=1.5)
    from_matrix = compute_harmonic_persistence(
        distance_matrix=distance_matrix(points, points), upper_bound=1.5
    )
    assert from_points["barcode"] == from_matrix["barcode"]
    assert from_points["harmonic"] == from_matrix["harmonic"]


def test_weights_sum_absolute_coefficients():
    points = two_circles()
    result = compute_harmonic_persistence(
        points=points, upper_bound=1.5, tolerance=0, vertex_labels=list("ab" * 15)
    )
    assert result["weights"]["vertices"] == list("ab" * 15)
    for representative, bar in zip(result["harmonic"], result["weights"]["bars"]):
        assert bar["bar_id"] == representative["bar_id"]
        expected = np.zeros(len(points))
        for term in representative["chain"]:
            expected[term["simplex"]] += abs(term["coefficient"])
        assert np.allclose(bar["weights"], expected)


def test_min_persistence_selects_long_bars():
    result = compute_harmonic_persistence(
        points=two_circles(), upper_bound=1.5, min_persistence=0.5
    )
    expected = [
        b["id"]
        for b in result["barcode"]
        if b["dimension"] == 1 and (b["persistence"] is None or b["persistence"] >= 0.5)
    ]
    assert expected
    assert [r["bar_id"] for r in result["harmonic"]] == expected


def test_barcode():
    result = compute_harmonic_persistence(points=two_circles(), upper_bound=1.5)
    # top-dimensional bars of a truncated Vietoris-Rips complex are left out
    assert {b["dimension"] for b in result["barcode"]} == {0, 1}
    for b in result["barcode"]:
        assert len(b["birth_simplex"]) == b["dimension"] + 1
        if b["death"] is None:
            assert b["persistence"] is None and b["death_simplex"] is None
        else:
            assert b["persistence"] == pytest.approx(b["death"] - b["birth"])
            assert len(b["death_simplex"]) == b["dimension"] + 2
    # the two circles are connected, and neither hole is filled below upper_bound
    infinite = [b["dimension"] for b in result["barcode"] if b["death"] is None]
    assert sorted(infinite) == [0, 1, 1]


def test_top_dimension_of_explicit_filtration():
    result = compute_harmonic_persistence(
        filtration=SPHERE_FILTRATION, homology_dimension=2
    )
    infinite = [b for b in result["barcode"] if b["death"] is None]
    assert sorted(b["dimension"] for b in infinite) == [0, 2]
    (representative,) = result["harmonic"]
    assert len(representative["chain"]) == 4
    assert len({abs(t["coefficient"]) for t in representative["chain"]}) == 1


def test_degree_zero():
    result = compute_harmonic_persistence(
        points=two_circles(), upper_bound=1.5, homology_dimension=0
    )
    assert len(result["harmonic"]) == len(two_circles())
    with pytest.raises(ValueError, match="precompute"):
        compute_harmonic_persistence(
            points=two_circles(),
            upper_bound=1.5,
            homology_dimension=0,
            precompute=False,
        )


def test_no_bars_in_selected_dimension():
    result = compute_harmonic_persistence(
        filtration=[([0], 0.0), ([1], 0.0), ([0, 1], 1.0)]
    )
    assert result["harmonic"] == [] and result["weights"]["bars"] == []


def test_filtration_order_only_matters_for_ties():
    # list order across different values and vertex order within simplices are
    # irrelevant; simplices with equal value and dimension keep their list order
    reordered = sorted(EXAMPLE_FILTRATION, key=lambda entry: -entry[1])
    reordered = [(s[::-1], v) for s, v in reordered]
    result = compute_harmonic_persistence(filtration=reordered)
    expected = compute_harmonic_persistence(filtration=EXAMPLE_FILTRATION)
    del result["timings"], expected["timings"]
    assert result == expected


def test_max_simplices():
    with pytest.raises(ValueError, match="more than the limit"):
        compute_harmonic_persistence(
            points=two_circles(), upper_bound=1.5, max_simplices=10
        )


@pytest.mark.parametrize(
    "kwargs, message",
    [
        (dict(), "exactly one"),
        (dict(points=[[0, 0]], distance_matrix=[[0]], upper_bound=1), "exactly one"),
        (dict(points=[[0, 0], [1, 1]]), "upper_bound' is required"),
        (dict(points=[[0, 0], [1, 1]], upper_bound=float("inf")), "positive finite"),
        (dict(points=[[0, 0], [1, np.nan]], upper_bound=1), "non-finite"),
        (dict(points=[0, 1], upper_bound=1), "2d array"),
        (dict(points=[[0, 0], [1, 1]], upper_bound=1, metric="nope"), "unknown metric"),
        (
            dict(points=[[0, 0], [1, 1]], upper_bound=1, homology_dimension=2),
            "smaller than 'max_dimension'",
        ),
        (dict(points=[[0, 0], [1, 1]], upper_bound=1, max_dimension=0), "at least 1"),
        (dict(distance_matrix=[[0, 1]], upper_bound=1), "square"),
        (dict(distance_matrix=[[0, 1], [2, 0]], upper_bound=1), "symmetric"),
        (dict(distance_matrix=[[1, 1], [1, 0]], upper_bound=1), "zero diagonal"),
        (dict(distance_matrix=[[0, -1], [-1, 0]], upper_bound=1), "negative"),
        (
            dict(points=[[0, 0], [1, 1]], upper_bound=1, vertex_labels=["a"]),
            "vertex labels",
        ),
        (dict(points=[[0, 0], [1, 1], [2, 2]], upper_bound=1, max_points=2), "limit"),
        (dict(filtration=[]), "empty"),
        (dict(filtration=[([0], 0), ([0, 1], 1)]), "face \\[1\\]"),
        (dict(filtration=[([0], 2), ([1], 0), ([0, 1], 1)]), "appears after"),
        (dict(filtration=[([0], 0), ([0], 1)]), "listed twice"),
        (dict(filtration=[([0], 0), ([2], 0)]), "0..n-1"),
        (dict(filtration=[([0, 0], 0)]), "repeats"),
        (dict(filtration=[([0], 0)], upper_bound=1), "only applies"),
        (dict(filtration=[([0], 0)], homology_dimension=1), "larger than"),
    ],
)
def test_invalid_inputs(kwargs, message):
    with pytest.raises(ValueError, match=message):
        validate_inputs(**kwargs)
