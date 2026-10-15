"""Harmonic persistent homology as a REST API.

:func:`matilda.api.service.compute_harmonic_persistence` runs the whole computation and
returns JSON-ready results; it only needs matilda's own dependencies. The FastAPI
application in :mod:`matilda.api.app` runs it as background jobs and needs the optional
dependencies, ``pip install "matilda[api]"``. Start the server with
``python -m matilda.api``.
"""

from .service import compute_harmonic_persistence, validate_inputs

__all__ = ["compute_harmonic_persistence", "validate_inputs"]
