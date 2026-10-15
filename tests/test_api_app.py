import time

import numpy as np
import pytest

pytest.importorskip("fastapi")
pytest.importorskip("httpx")

from fastapi.testclient import TestClient

from matilda.api import compute_harmonic_persistence
from matilda.api.app import Settings, create_app

POINTS = [[np.cos(t), np.sin(t)] for t in np.linspace(0, 2 * np.pi, 12, endpoint=False)]
REQUEST = {"points": POINTS, "upper_bound": 1.2, "vertex_labels": list("abcdefghijkl")}


@pytest.fixture
def make_client():
    clients = []

    def make(**settings):
        client = TestClient(create_app(Settings(**settings)))
        client.__enter__()
        clients.append(client)
        return client

    yield make
    for client in clients:
        client.__exit__(None, None, None)


def wait(client, job_id, timeout=120):
    deadline = time.time() + timeout
    while time.time() < deadline:
        status = client.get("/jobs/{}".format(job_id)).json()
        if status["status"] in ("succeeded", "failed", "cancelled"):
            return status
        time.sleep(0.1)
    raise AssertionError("job {} did not finish".format(job_id))


def test_job_lifecycle(make_client):
    client = make_client()
    assert client.get("/health").json() == {"status": "ok"}

    response = client.post("/jobs", json=REQUEST)
    assert response.status_code == 202
    job = response.json()
    assert job["status"] in ("queued", "running")

    status = wait(client, job["id"])
    assert status["status"] == "succeeded", status["error"]
    assert status["summary"]["n_simplices"][0] == len(POINTS)
    assert status["summary"]["n_representatives"] == 1

    expected = compute_harmonic_persistence(**REQUEST)
    barcode = client.get("/jobs/{}/barcode".format(job["id"])).json()
    assert barcode["bars"] == expected["barcode"]
    harmonic = client.get("/jobs/{}/harmonic".format(job["id"])).json()
    assert harmonic["homology_dimension"] == 1
    assert harmonic["representatives"] == expected["harmonic"]
    weights = client.get("/jobs/{}/weights".format(job["id"])).json()
    assert weights["vertices"] == list("abcdefghijkl")
    assert weights["bars"] == expected["weights"]["bars"]

    assert [j["id"] for j in client.get("/jobs").json()] == [job["id"]]
    assert client.delete("/jobs/{}".format(job["id"])).status_code == 204
    assert client.get("/jobs/{}".format(job["id"])).status_code == 404


def test_filtration_job(make_client):
    client = make_client()
    request = {
        "filtration": [
            {"simplex": s, "value": v}
            for s, v in [
                ([0], 0),
                ([1], 0),
                ([2], 0),
                ([0, 1], 1),
                ([1, 2], 1),
                ([0, 2], 2),
            ]
        ]
    }
    job = client.post("/jobs", json=request).json()
    assert wait(client, job["id"])["status"] == "succeeded"
    (representative,) = client.get("/jobs/{}/harmonic".format(job["id"])).json()[
        "representatives"
    ]
    assert representative["death"] is None
    assert sorted(t["simplex"] for t in representative["chain"]) == [
        [0, 1],
        [0, 2],
        [1, 2],
    ]


@pytest.mark.parametrize(
    "request_json, message",
    [
        ({"points": POINTS}, "upper_bound"),
        ({"points": POINTS, "upper_bound": 1, "distance_matrix": [[0]]}, "exactly one"),
        ({"points": POINTS, "upper_bound": 1, "unknown": 1}, "unknown"),
        ({"points": POINTS, "upper_bound": -1}, "upper_bound"),
        ({"filtration": [{"simplex": [0, 1], "value": 0}]}, "missing"),
    ],
)
def test_invalid_request(make_client, request_json, message):
    response = make_client().post("/jobs", json=request_json)
    assert response.status_code == 422
    assert message in str(response.json()["detail"])


def test_max_points(make_client):
    response = make_client(max_points=5).post("/jobs", json=REQUEST)
    assert response.status_code == 422
    assert "limit of 5" in response.json()["detail"]


def test_unknown_job(make_client):
    client = make_client()
    for path in ("/jobs/nope", "/jobs/nope/barcode", "/jobs/nope/weights"):
        assert client.get(path).status_code == 404
    assert client.delete("/jobs/nope").status_code == 404


def test_failed_job(make_client):
    client = make_client(max_simplices=10)
    job = client.post("/jobs", json=REQUEST).json()
    status = wait(client, job["id"])
    assert status["status"] == "failed"
    assert "more than the limit of 10" in status["error"]
    response = client.get("/jobs/{}/harmonic".format(job["id"]))
    assert response.status_code == 409
    assert "more than the limit of 10" in response.json()["detail"]


def test_job_timeout(make_client):
    client = make_client(job_timeout=0.01)
    job = client.post("/jobs", json=REQUEST).json()
    status = wait(client, job["id"])
    assert status["status"] == "failed"
    assert "time limit" in status["error"]


def test_deleting_a_running_job_frees_its_worker(make_client):
    client = make_client(max_workers=1)
    first = client.post("/jobs", json=REQUEST).json()
    assert client.delete("/jobs/{}".format(first["id"])).status_code == 204
    second = client.post("/jobs", json=REQUEST).json()
    assert wait(client, second["id"])["status"] == "succeeded"
