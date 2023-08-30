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