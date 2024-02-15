import numpy as np
import pandas as pd
import networkx as nx

import time
import matplotlib.pyplot as plt

import matilda
import matilda.prototyping
import matilda.harmonic

from scipy.stats import median_abs_deviation


def MAD_normalization(X):
    # The expression data are also normalized using the median absolute deviation (MAD) scale normalization,
    # which is calculated as (x-median)/MAD, where x is the gene expression data.
    # https://omicssimla.sourceforge.io/RNASeq_norgene.html

    return (X - np.median(X, axis=0)) / median_abs_deviation(X)


def cycles_longer_that(homology_computer, min_len, upper_bound, H=1):
    """returns cycles whose bars are longer than `min_len`

    Parameters
    ----------
    homology_computer : matilda.PersistentHomologyComputer()
        fitted matilda homology computer
    min_len : float
        the min bar lenght
    upper_bound : float
        the upper bound used in `construct_vietoris_from_metric`
    H : int, optional
        cycle dimension, by default 1

    Returns
    -------
    dict
        dictionary with the cycles longer than `min_len`
    """

    return {
        k: c
        for k, c in homology_computer.persistent_cycles[H].items()
        if min(homology_computer.bars[H][k][1], upper_bound)
        - homology_computer.bars[H][k][0]
        > min_len
    }


def sorted_cycles(homology_computer, upper_bound, H=1):
    cycles_lenght = {}

    for k in homology_computer.persistent_cycles[H]:
        cycles_lenght[k] = (
            min(homology_computer.bars[H][k][1], upper_bound)
            - homology_computer.bars[H][k][0]
        )

    return {
        k: v
        for k, v in sorted(
            cycles_lenght.items(), key=lambda item: item[1], reverse=True
        )
    }


def get_harmonic_weights(my_harmonic_cycles, K, index_list):
    """takes as input a dict of harmonic cycles and returns a dict of weights on the simplices and on the nodes

    Parameters
    ----------
    my_harmonic_cycles : dict
        dict of harmonic cycles, as per matilda output
    K : matilda.FilteredSimplicialComplex
        the simplicial complex whose homology was computed
    index_list : list
        list of node names (the samples ids)

    Returns
    -------
    (dict, dict)
        dict of weights on the nodes, dict of weights on the simplices
    """

    simplices_dict = {}
    nodes_dict = {}

    for k, h in my_harmonic_cycles.items():
        simplices_dict[k] = {
            i: {"simplex": K.simplices[K.simplices_indices[i]], "weight": np.abs(w)}
            for i, w in h.items()
        }

        nodes_dict[k] = {}

        for i, s in simplices_dict[k].items():
            for n in s["simplex"]:
                nodes_dict[k][n] = nodes_dict[k].get(n, dict())

                nodes_dict[k][n]["weight"] = (
                    nodes_dict[k][n].get("weight", 0) + s["weight"]
                )
                nodes_dict[k][n]["id"] = index_list[n]

    return nodes_dict, simplices_dict


def get_harmonic_weights_df(
    nodes_dict, my_cycles, index_list, homology_computer, upper_bound, H=1
):
    weights_df = pd.DataFrame(
        np.array(
            [
                [
                    nodes_dict[k].get(i, {"weight": 0})["weight"]
                    for i, _ in enumerate(index_list)
                ]
                for k in nodes_dict
            ]
        ),
        columns=index_list,
    )

    weights_df["bar_id"] = [i for i in my_cycles]
    weights_df["birth"] = [homology_computer.bars[H][i][0] for i in my_cycles]
    weights_df["death"] = [homology_computer.bars[H][i][1] for i in my_cycles]
    weights_df["bar_length"] = [
        min(homology_computer.bars[H][i][1], upper_bound)
        - homology_computer.bars[H][i][0]
        for i in my_cycles
    ]

    weights_df = weights_df[
        weights_df.columns[-4:].tolist() + weights_df.columns[:-4].tolist()
    ]

    return weights_df


def pipeline(
    dm,
    max_dimension,
    upper_bound,
    cycle_dimension,
    gene_names,
    show_figures=False,
    **kwargs
):
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

    if show_figures:
        plotter = matilda.plot.Plotter()
        fig, ax = plotter.plot_barcode(
            homology_computer,
            dimension=[cycle_dimension],
            figsize=(4, 4),
        )
        fig.tight_layout()
        plt.show()
    print("{:.2f}s".format(time.time() - start))

    harmonic_computer = matilda.harmonic.HarmonicRepresentativesComputer(
        K, homology_computer
    )

    start = time.time()
    print(
        "projecting {} {}-dimensional cycles".format(
            len(homology_computer.bars[cycle_dimension]), cycle_dimension
        )
    )

    harmonic_computer.compute_harmonic_cycles(
        dim=cycle_dimension, precompute=True, verbose=1
    )
    print("{:.2f}s".format(time.time() - start))

    # return {
    #     "bars": homology_computer.bars,
    #     "cycles": homology_computer.persistent_cycles,
    #     "harmonic_cycles": harmonic_computer.harmonic_cycles,
    # }

    # let's record all weights of simplices in any harmonic repr
    normalized_weights_dict = {}

    print("there are {} 1-d cycles".format(len(harmonic_computer.harmonic_cycles[1])))
    total_bar_len = 0
    for c in harmonic_computer.harmonic_cycles[1]:
        total_bar_len += (
            homology_computer.bars[1][c][1] - homology_computer.bars[1][c][0]
        )

    for c in harmonic_computer.harmonic_cycles[1]:
        len_bar = homology_computer.bars[1][c][1] - homology_computer.bars[1][c][0]

        for k, w in harmonic_computer.harmonic_cycles[1][c].items():
            normalized_weights_dict[k] = (
                normalized_weights_dict.get(k, 0) + np.abs(w) * len_bar / total_bar_len
            )

    normalized_weights_dict = dict(
        sorted(normalized_weights_dict.items(), key=lambda item: item[1], reverse=True)
    )

    edges_list = [
        (*K.simplices[i], {"weight": w}) for i, w in normalized_weights_dict.items()
    ]

    G = nx.Graph()
    G.add_edges_from(edges_list)

    del harmonic_computer
    del K

    for n in G.nodes:
        G.nodes[n]["weight"] = sum([e[2]["weight"] for e in G.edges(n, data=True)])
        G.nodes[n]["gene"] = gene_names[n]

    # aggreggate weights on the nodes
    return {n[1]["gene"]: n[1]["weight"] for n in G.nodes(data=True)}


def load_data(filenames, include_DE=True):
    df_list = []
    for i, f in enumerate(filenames):
        file_type = f.split(".")[-1]
        if i > 0:
            assert file_type == first_file_type
        else:
            first_file_type = file_type
        # rna
        if file_type == "exp":
            df = pd.read_csv(f, delimiter=" ").dropna(axis=1)
            df["is_case"] = [int("DCASE" in i) for i in df.index]
            df = df[["is_case"] + df.columns[:-1].tolist()]

        # protein
        elif file_type == "pro":
            df = pd.read_csv(f, delimiter=" ", header=None, index_col=0).dropna(axis=1)
            df.index.name = None
            df["is_case"] = [int("DCASE" in i) for i in df.index]
            df.rename(
                {
                    i: n
                    for i, n in enumerate(
                        ["FAM", "ID", "F_ID", "M_ID", "SEX", "AFF"], start=1
                    )
                },
                axis=1,
                inplace=True,
            )
            df = df[["is_case"] + df.columns[:-1].tolist()]

        # methylation
        elif file_type == "methy":
            df = pd.read_csv(f, delimiter=" ", header=None, index_col=0).dropna(axis=1)
            df.index.name = None
            df["is_case"] = [int("DCASE" in i) for i in df.index]
            df = df[["is_case"] + df.columns[:-1].tolist()]

        df_list.append(df)

    full_df = pd.concat(df_list)

    if file_type == "methy":
        X = full_df[full_df.columns[3:]]
        X = X.loc[:, (X != "-1,-1").all(axis=0)]
        X = X.map(lambda x: int(x.split(",")[0]) / int(x.split(",")[1]))

    elif include_DE or file_type != "exp":
        X = full_df[full_df.columns[7:]]

    else:
        X = full_df[full_df.columns[8:]]

    y = full_df.is_case.to_numpy(dtype=int)

    return X, y


def plot_bars(genes_weights, genes_to_plot):
    fig, ax = plt.subplots(figsize=(10, 4), ncols=1, nrows=1)

    labels_dict = {
        0: "control",
        1: "case",
    }

    x = np.arange(len(genes_to_plot))  # the label locations
    width = 0.25  # the width of the bars
    multiplier = 0

    for t in genes_weights:
        offset = width * multiplier
        rects = ax.bar(
            x + offset,
            [genes_weights[t][e] for e in genes_to_plot],
            width,
            label=labels_dict[t],
        )
        multiplier += 1

        ax.set_xticks(x + width, genes_to_plot, rotation=45)
        # ax.xaxis.set_tick_params(length=20)

    plt.legend()
    plt.tight_layout()
    # plt.savefig('brca_SCIPY.pdf')
    return fig, ax
