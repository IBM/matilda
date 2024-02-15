import numpy as np
import pandas as pd
import scipy

from tqdm.auto import tqdm

import matplotlib.pyplot as plt
import networkx as nx

import matilda
import matilda.prototyping
import matilda.harmonic
import seaborn as sns

import time

from sklearn.metrics import pairwise_distances

import utils


def extract_harmonic_weights(
    X,
    upper_bound,
    metric="euclidean",
    output_file=None,
    max_dimension=2,
    H=1,
    show_figures=True,
):
    """pipeline to extract the harmonic weights

    Parameters
    ----------
    X : numpy.array
        the dataset
    upper_bound : float
        max filtration radius for the VR complex
    metric : string, optional
        metric to use to compute the distance matrix, by default 'euclidean'
    output_file : string, optional
        path of the pickle file where to save the output, by default None
    max_dimension : int, optional
        max dimension of the VR complex, by default 2
    H : int, optional
        dimension of the homology to compute, by default 1
    show_figures : bool, optional
        whether to display figures, by default True

    Returns
    -------
    pandas.DataFrame
        a DataFrame with the harmonic weights for each bar for each input datapoint
    """
    ## COMPUTE DISTANCE MATRIX
    print("computing distance matrix")

    dm = pairwise_distances(X=X, metric=metric, n_jobs=-1)
    np.fill_diagonal(dm, 0)

    print(dm.shape)

    if show_figures:
        fig, axs = plt.subplots(ncols=2, figsize=(10, 4))
        axs[0].imshow(dm)
        sns.histplot(dm.ravel(), ax=axs[1])
        plt.suptitle(metric)
        plt.show()

    ## HOMOLOGY
    start = time.time()
    # create simplicial complex
    K_frompoints_co = matilda.FilteredSimplicialComplex()
    K_frompoints_co.construct_vietoris_from_metric(
        dm, dimension=max_dimension, upper_bound=upper_bound
    )

    # we need to use the prototype class
    K = matilda.prototyping.FilteredSimplicialComplex(
        dimension=K_frompoints_co.dimension,
        simplices=K_frompoints_co.simplices,
        simplices_indices=K_frompoints_co.simplices_indices,
        appears_at=K_frompoints_co.appears_at,
    )

    print("complex created {:.2f}s".format(time.time() - start))
    for i in range(K.dimension + 1):
        print(
            "{} {}-dim simplices".format(
                len([s for s in K.simplices if len(s) == (i + 1)]), i
            )
        )

    print("computing homology")
    start = time.time()
    homology_computer = matilda.PersistentHomologyComputer()
    homology_computer.compute_persistent_homology(
        K, with_representatives=True, modulus=0
    )

    print("done {:.2f}s".format(time.time() - start))

    if show_figures:
        plotter = matilda.plot.Plotter()
        fig, ax = plotter.plot_barcode(
            homology_computer, dimension=[H], figsize=(4, 4), max_x=upper_bound + 0.05
        )
        fig.tight_layout()
        plt.show()

    print("there are {} {}-dimensional bars".format(len(homology_computer.bars[H]), H))

    ## HARMONICS
    print("computing harmonics")
    start = time.time()

    harmonic_computer = matilda.harmonic.HarmonicRepresentativesComputer(
        K, homology_computer
    )

    harmonic_computer.compute_harmonic_cycles(dim=H, verbose=1, precompute=True)
    print("done {:.2f}s".format(time.time() - start))

    nodes_dict, simplices_dict = utils.get_harmonic_weights(
        harmonic_computer.harmonic_cycles[H], K, index_list=X.index.to_list()
    )

    weights_df = utils.get_harmonic_weights_df(
        nodes_dict,
        homology_computer.persistent_cycles[H],
        X.index.to_list(),
        homology_computer,
        upper_bound,
    )

    weights_df.set_index("bar_id", inplace=True)
    weights_df.sort_values("bar_length", inplace=True, ascending=False)

    ## SAVING
    if output_file:
        weights_df.to_pickle(output_file)

    return weights_df
