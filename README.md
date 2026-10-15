# maTilDA - Multipurpose toolkit for TDA
## Introduction
**maTilDA** is a library aimed towards topological data analysis researchers and practicioners with emphasis on available features, speed, and ease of use.

Currently, **maTilDA** supports creation, and some manipulations, of filtered simplicial complexes, along with persistent homology computations over $\mathbb{F}_p$ or $\mathbb{R}$. Filtered complexes can be constructed from euclidean spaces or more generally from metric spaces. Homology computations include usual barcode/persistence diagram representations, recovery of representative cycles, and derived summaries.

Additionally, **maTilDA** also includes mapper functionality, with support of standard rectangular coverings and custom coverings. Mapper outputs can be visualized in an interactive manner through a standalone html file or as output of a jupyter notebook cell.

## Installation

1. Download this repository, either cloning with git or using the "Download Zip" option available on the repository's github webpage.
2. Assuming you downloaded the zip file, you can install matilda by running 

        pip install ./matilda-master.zip
    from the folder containing the zip file.
3. That's it! All dependencies should be installed automatically.

## Quickstart

We will walk you through the notebook `example_persistence.ipynb`. 

The first cell starts with:

```python
import matilda
``` 
which, naturally, imports the module. This is followed by

```python
matrix = matilda.examples.circle_metric(res=100)
```
**maTilDA** includes several ready-to-use toy examples presented as distance matrices, or as point clouds in euclidean space. In this case we instance a distance matrix of 100 points lying in a unitary circle and assign it to `matrix`. Next we have

```python

K = matilda.FilteredSimplicialComplex()
K.construct_vietoris_from_metric(matrix,2,2.1)
```

The class `FilteredSimplicialComplex` is the main class used for persistent homology computations. In the code above an empty filtered simplicial complex was created, and then, by using the `construct_vietoris_from_metric` method, we fill it up with a Vietoris-Rips complex, constructed from the data in `matrix`. In this case, we construct all simplices up to dimension `2` and we consider only pairs of points with a distance smaller than `2.1`. 

After this, we create an instance of `PersistentHomologyComputer`, which is the class used to perform homology computations. In particular it has the method `compute_persistent_homology`, which, perhaps unsurprisingly, receives as a parameter a filtered simplicial complex and computes its persistent homology.

```python
homology_computer = matilda.PersistentHomologyComputer()

homology_computer.compute_persistent_homology(K)
```

After the computation is finished, the object `homology_computer` now has the persistent homology information of `K`. We can inspect it manually, the way we would inspect a python object, or we can produce a visualization. To do this, we create an instance of `Plotter`, a class to produce persistent homology visualizations. Here we use the `plot_barcode` method, which takes as input a homology computer object and optionally the homology degrees of interest.

```python
plotter = matilda.plot.Plotter()
plotter.plot_barcode(homology_computer, dimension=[0,1])
```

The last cell contains a way to print out the bars, which are stored as a python dictionary in the `bars` attribute of the plotter object. Note that since this is a complex of dimension `2`, there are no boundaries in degree `2` and therefore all cycles will be present as bars. To prevent them from being presented and filling up the cell's output, we print only bars of degrees `0` and `1`.

```python
for k,v in homology_computer.bars.items():
    if k>1:
        break
    print("Bars of dimension {}".format(k))
    for kk,vv in v.items():
        print("{}:{}".format(kk,vv))
```

## REST API for harmonic persistent homology

**maTilDA** can serve harmonic persistent homology over HTTP. The API computes the barcode of a filtered complex over $\mathbb{R}$ and the harmonic representatives of its bars with `matilda.harmonic.HarmonicRepresentativesComputer`, together with the per-vertex harmonic weights used in Gurnari et al., "Probing omics data via harmonic persistent homology".

### Running the server

Install with the optional API dependencies and start the server:

    pip install "./matilda-master.zip[api]"
    python -m matilda.api --host 127.0.0.1 --port 8000

Interactive documentation is then available at `http://127.0.0.1:8000/docs`. The notebook `examples/example_api.ipynb` walks through the API, computing barcodes and harmonic representatives over HTTP. The server has no authentication, so it binds to `127.0.0.1` by default; do not expose it on a public network.

### Submitting a job

Computations can take minutes, so each request runs as a background job in its own process. Submit a job, poll its status until it is `succeeded`, then fetch the results:

```python
import time
import numpy as np
import requests

points = np.random.default_rng(0).normal(size=(100, 2))
job = requests.post(
    "http://127.0.0.1:8000/jobs",
    json={"points": points.tolist(), "upper_bound": 1.0, "homology_dimension": 1},
).json()

while job["status"] in ("queued", "running"):
    time.sleep(1)
    job = requests.get("http://127.0.0.1:8000/jobs/" + job["id"]).json()

bars = requests.get("http://127.0.0.1:8000/jobs/{}/barcode".format(job["id"])).json()["bars"]
harmonic = requests.get("http://127.0.0.1:8000/jobs/{}/harmonic".format(job["id"])).json()
weights = requests.get("http://127.0.0.1:8000/jobs/{}/weights".format(job["id"])).json()
```

| Endpoint | Description |
|---|---|
| `POST /jobs` | validate the input and queue a job; returns the job status with its `id` |
| `GET /jobs` | status of all jobs kept in memory |
| `GET /jobs/{id}` | status of a job: `queued`, `running`, `succeeded`, `failed` (with an `error`) or `cancelled` |
| `GET /jobs/{id}/barcode` | all bars, each with an `id`, `dimension`, `birth`, `death`, `persistence` and the simplices creating and killing it |
| `GET /jobs/{id}/harmonic` | harmonic representatives of the selected bars, as chains of `{"simplex", "coefficient"}` terms |
| `GET /jobs/{id}/weights` | for each selected bar, the weight of every vertex: the sum of the absolute coefficients of the simplices containing it |
| `DELETE /jobs/{id}` | cancel a job if it is still running and forget it |

The body of `POST /jobs` gives exactly one input complex:

- `points`: an $n \times d$ point cloud, with `metric` (default `"euclidean"`) and `upper_bound`, the maximum edge length of the Vietoris-Rips complex;
- `distance_matrix`: an $n \times n$ distance matrix, with `upper_bound`;
- `filtration`: an explicit filtered complex, as a list of `{"simplex": [vertices], "value": filtration value}`. Every face must be listed, with a value not larger than its cofaces', and vertices are numbered `0..n-1`. Simplices with equal values are ordered by dimension and then by their position in the list.

Other fields are `max_dimension` (maximum simplex dimension of the Vietoris-Rips complex, default 2), `homology_dimension` (default 1), `min_persistence` (only bars at least this long get representatives; infinite bars always do), `precompute` (passed to `compute_harmonic_cycles`, default `true`), `tolerance` (chain coefficients at most this in absolute value are dropped, default `1e-12`) and `vertex_labels` (one label per vertex, used in the weights output). For Vietoris-Rips inputs the barcode leaves out the top dimension, where every cycle is an infinite bar.

Bar ids are opaque and only meaningful within a job. Simplices are given by their vertices.

### Limits

Jobs and results are kept in memory. The server reads its limits from environment variables; set a variable to `none` to remove a limit.

| Variable | Default | Description |
|---|---|---|
| `MATILDA_API_MAX_WORKERS` | 1 | jobs running concurrently |
| `MATILDA_API_JOB_TIMEOUT` | 3600 | seconds after which a running job is killed |
| `MATILDA_API_MAX_MEMORY_MB` | none | memory limit of each job's process |
| `MATILDA_API_MAX_POINTS` | 5000 | maximum number of vertices of an input |
| `MATILDA_API_MAX_SIMPLICES` | 5000000 | maximum number of simplices of a complex |
| `MATILDA_API_MAX_FINISHED_JOBS` | 100 | finished jobs kept; the oldest are forgotten first |

### Using it from Python

The same computation is available without the server, and without the API dependencies:

```python
from matilda.api import compute_harmonic_persistence

result = compute_harmonic_persistence(points=points, upper_bound=1.0)
result["barcode"], result["harmonic"], result["weights"]
```

The API's tests run with `pip install "./matilda-master.zip[api,test]"` and `pytest tests`.
