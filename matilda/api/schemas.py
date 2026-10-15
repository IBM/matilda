"""Request and response models of the harmonic persistence API."""

from typing import Dict, List, Optional, Union

from pydantic import BaseModel, ConfigDict, Field


class FiltrationEntry(BaseModel):
    simplex: List[int] = Field(min_length=1, description="vertices of the simplex")
    value: float = Field(description="filtration value at which the simplex appears")


class JobRequest(BaseModel):
    """input complex and parameters of a harmonic persistence job

    Give exactly one of `points`, `distance_matrix` or `filtration`.
    """

    model_config = ConfigDict(extra="forbid")

    points: Optional[List[List[float]]] = Field(
        None, description="(n, d) point cloud; builds a Vietoris-Rips complex"
    )
    distance_matrix: Optional[List[List[float]]] = Field(
        None, description="(n, n) distance matrix; builds a Vietoris-Rips complex"
    )
    filtration: Optional[List[FiltrationEntry]] = Field(
        None,
        description="explicit filtered complex; every face must be listed, "
        "vertices numbered 0..n-1. Simplices are ordered by value, then dimension, "
        "then position in this list",
    )
    metric: str = Field("euclidean", description="distance used for `points`")
    max_dimension: int = Field(
        2, ge=1, description="maximum simplex dimension of the Vietoris-Rips complex"
    )
    upper_bound: Optional[float] = Field(
        None,
        gt=0,
        description="maximum edge length of the Vietoris-Rips complex; "
        "required for `points` and `distance_matrix`",
    )
    homology_dimension: int = Field(
        1,
        ge=0,
        description="dimension of the bars to compute harmonic representatives of",
    )
    min_persistence: float = Field(
        0.0,
        ge=0,
        description="only return representatives of bars at least this long; "
        "infinite bars are always returned",
    )
    precompute: bool = Field(
        True, description="orthonormalize the whole boundary matrix once"
    )
    tolerance: float = Field(
        1e-12,
        ge=0,
        description="drop chain coefficients at most this in absolute value",
    )
    vertex_labels: Optional[List[str]] = Field(
        None, description="one label per vertex, used in the weights output"
    )

    def to_params(self):
        """keyword arguments of `matilda.api.service.compute_harmonic_persistence`"""
        params = self.model_dump()
        if self.filtration is not None:
            params["filtration"] = [(e.simplex, e.value) for e in self.filtration]
        return params


class JobSummary(BaseModel):
    input: Dict[str, Union[str, int, float, bool, None]]
    n_simplices: List[int]
    n_bars: int
    n_representatives: int
    timings: Dict[str, float]


class JobStatus(BaseModel):
    id: str
    status: str = Field(description="queued, running, succeeded, failed or cancelled")
    created_at: str
    started_at: Optional[str]
    finished_at: Optional[str]
    error: Optional[str]
    summary: Optional[JobSummary] = Field(
        description="present once the job has succeeded"
    )


class Bar(BaseModel):
    id: int = Field(description="identifier of the bar, unique within a job")
    dimension: int
    birth: float
    death: Optional[float] = Field(description="None for infinite bars")
    persistence: Optional[float] = Field(description="None for infinite bars")
    birth_simplex: List[int] = Field(
        description="simplex whose addition creates the bar"
    )
    death_simplex: Optional[List[int]] = Field(
        description="simplex whose addition kills the bar, None for infinite bars"
    )


class Barcode(BaseModel):
    bars: List[Bar]


class ChainTerm(BaseModel):
    simplex: List[int]
    coefficient: float


class HarmonicRepresentative(BaseModel):
    bar_id: int
    birth: float
    death: Optional[float]
    persistence: Optional[float]
    chain: List[ChainTerm]


class HarmonicRepresentatives(BaseModel):
    homology_dimension: int
    vertex_labels: Optional[List[str]]
    representatives: List[HarmonicRepresentative]


class BarWeights(BaseModel):
    bar_id: int
    birth: float
    death: Optional[float]
    persistence: Optional[float]
    weights: List[float] = Field(
        description="one weight per vertex, in `vertices` order"
    )


class VertexWeights(BaseModel):
    homology_dimension: int
    vertices: List[Union[str, int]] = Field(
        description="vertex labels if given, else vertex indices"
    )
    bars: List[BarWeights]
