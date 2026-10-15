"""FastAPI application exposing harmonic persistent homology as background jobs.

Submit a job with ``POST /jobs``, poll ``GET /jobs/{id}`` until its status is
"succeeded", then fetch ``/jobs/{id}/barcode``, ``/jobs/{id}/harmonic`` or
``/jobs/{id}/weights``.
"""

import os
from contextlib import asynccontextmanager
from dataclasses import dataclass
from typing import List, Optional

from fastapi import FastAPI, HTTPException, Response
from fastapi.responses import JSONResponse

from .jobs import SUCCEEDED, JobManager
from .schemas import (
    Barcode,
    HarmonicRepresentatives,
    JobRequest,
    JobStatus,
    VertexWeights,
)
from .service import validate_inputs


def _env(name, convert, default):
    value = os.environ.get(name)
    if value is None or value == "":
        return default
    if value.lower() == "none":
        return None
    return convert(value)


@dataclass
class Settings:
    """limits of the API; `from_env` reads them from MATILDA_API_* variables

    Attributes
    ----------
    max_workers : int
        jobs running concurrently (MATILDA_API_MAX_WORKERS)
    job_timeout : float or None
        seconds after which a running job is killed (MATILDA_API_JOB_TIMEOUT)
    max_memory_mb : int or None
        address-space limit of each job's subprocess (MATILDA_API_MAX_MEMORY_MB)
    max_points : int or None
        maximum number of vertices of an input (MATILDA_API_MAX_POINTS)
    max_simplices : int or None
        maximum number of simplices of a complex (MATILDA_API_MAX_SIMPLICES)
    max_finished_jobs : int
        finished jobs kept in memory (MATILDA_API_MAX_FINISHED_JOBS)
    """

    max_workers: int = 1
    job_timeout: Optional[float] = 3600.0
    max_memory_mb: Optional[int] = None
    max_points: Optional[int] = 5000
    max_simplices: Optional[int] = 5_000_000
    max_finished_jobs: int = 100

    @classmethod
    def from_env(cls):
        default = cls()
        return cls(
            max_workers=_env("MATILDA_API_MAX_WORKERS", int, default.max_workers),
            job_timeout=_env("MATILDA_API_JOB_TIMEOUT", float, default.job_timeout),
            max_memory_mb=_env("MATILDA_API_MAX_MEMORY_MB", int, default.max_memory_mb),
            max_points=_env("MATILDA_API_MAX_POINTS", int, default.max_points),
            max_simplices=_env("MATILDA_API_MAX_SIMPLICES", int, default.max_simplices),
            max_finished_jobs=_env(
                "MATILDA_API_MAX_FINISHED_JOBS", int, default.max_finished_jobs
            ),
        )


def create_app(settings=None):
    """build the FastAPI application

    Parameters
    ----------
    settings : Settings, optional
        limits of the API, by default read from the environment with `Settings.from_env`

    Returns
    -------
    fastapi.FastAPI
        the application; its job manager is `app.state.jobs`
    """
    if settings is None:
        settings = Settings.from_env()
    jobs = JobManager(
        max_workers=settings.max_workers,
        job_timeout=settings.job_timeout,
        max_memory_mb=settings.max_memory_mb,
        max_finished_jobs=settings.max_finished_jobs,
    )

    @asynccontextmanager
    async def lifespan(app):
        yield
        jobs.shutdown()

    app = FastAPI(
        title="matilda harmonic persistent homology",
        description="Barcodes and harmonic representatives of filtered simplicial "
        "complexes, computed with matilda as background jobs.",
        lifespan=lifespan,
    )
    app.state.settings = settings
    app.state.jobs = jobs

    def get_job(job_id):
        job = jobs.get(job_id)
        if job is None:
            raise HTTPException(status_code=404, detail="unknown job {}".format(job_id))
        return job

    def get_result(job_id):
        job = get_job(job_id)
        if job.status != SUCCEEDED:
            detail = "job {} is {}".format(job_id, job.status)
            if job.error:
                detail += ": " + job.error
            raise HTTPException(status_code=409, detail=detail)
        return job.result

    @app.get("/health")
    def health():
        return {"status": "ok"}

    @app.post("/jobs", status_code=202, response_model=JobStatus)
    def submit_job(request: JobRequest):
        """validate the input and queue a job computing harmonic persistence"""
        params = request.to_params()
        try:
            checked = validate_inputs(
                **{
                    name: params[name]
                    for name in (
                        "points",
                        "distance_matrix",
                        "filtration",
                        "metric",
                        "max_dimension",
                        "upper_bound",
                        "homology_dimension",
                        "vertex_labels",
                    )
                },
                max_points=settings.max_points,
            )
        except ValueError as e:
            raise HTTPException(status_code=422, detail=str(e))
        # numpy arrays are much cheaper than nested lists to send to the job's process
        params.update(checked, max_simplices=settings.max_simplices)
        return jobs.submit(params).describe()

    @app.get("/jobs", response_model=List[JobStatus])
    def list_jobs():
        """all jobs kept in memory, oldest first"""
        return [job.describe() for job in jobs.list_jobs()]

    @app.get("/jobs/{job_id}", response_model=JobStatus)
    def job_status(job_id: str):
        return get_job(job_id).describe()

    @app.delete("/jobs/{job_id}", status_code=204)
    def delete_job(job_id: str):
        """cancel the job if it is still running and forget it"""
        if not jobs.delete(job_id):
            raise HTTPException(status_code=404, detail="unknown job {}".format(job_id))
        return Response(status_code=204)

    # results can hold millions of numbers, so they are returned as JSONResponse to skip
    # re-validating them against the response models, which only document the shape

    @app.get("/jobs/{job_id}/barcode", response_model=Barcode)
    def barcode(job_id: str):
        """all bars, in every dimension up to the meaningful one"""
        return JSONResponse({"bars": get_result(job_id)["barcode"]})

    @app.get("/jobs/{job_id}/harmonic", response_model=HarmonicRepresentatives)
    def harmonic(job_id: str):
        """harmonic representatives of the selected bars, as chains of simplices"""
        result = get_result(job_id)
        return JSONResponse(
            {
                "homology_dimension": result["input"]["homology_dimension"],
                "vertex_labels": result["vertex_labels"],
                "representatives": result["harmonic"],
            }
        )

    @app.get("/jobs/{job_id}/weights", response_model=VertexWeights)
    def weights(job_id: str):
        """per-vertex harmonic weights of the selected bars"""
        result = get_result(job_id)
        return JSONResponse(
            {
                "homology_dimension": result["input"]["homology_dimension"],
                **result["weights"],
            }
        )

    return app
