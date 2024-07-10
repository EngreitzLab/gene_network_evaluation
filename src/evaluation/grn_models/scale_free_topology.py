import os
import argparse

import mudata

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from igraph import Graph
import networkx as nx
from sklearn.linear_model import LinearRegression


def _compute_scale_free_topology_fit_deprecated(
    grn: pd.DataFrame,
    source: str = "tf",
    target: str = "gene",
    score: str = "score",
):

    # Build iGraph
    g = Graph.DataFrame(grn[[source, target]], directed=True, use_vids=False)
    g.es[score] = grn[score].values.copy()

    # Extract degree distribution
    degree_df = pd.DataFrame(g.degree(mode="all"), columns=["degree"])
    dist = degree_df.degree.value_counts()/degree_df.degree.value_counts().sum()
    
    # Fit power law
    x = np.log(dist.index.values).reshape([-1,1])
    y = np.log(dist.values).reshape([-1,1])
    model = LinearRegression()
    model.fit(x,y)

    # Store degree distribution, model, and R2 in uns["scale_free_topology"]
    output = {
        "degree_distribution": dist,
        "model": model,
        "r2": model.score(x,y)
    }

    return output


def _compute_scale_free_topology_fit(
    grn: pd.DataFrame,
    source: str = "tf",
    target: str = "gene",
):
    # Make networkx graph
    G = nx.DiGraph()
    G_ = nx.from_pandas_edgelist(grn, source=source, target=target, edge_attr=True)
    G.add_edges_from(G_.edges())

    # Get degree distribution
    df = pd.DataFrame(np.array(G.degree))
    df.columns = ["ind", "degree"]
    df = df.set_index("ind")
    dist = df.degree.value_counts()/df.degree.value_counts().sum()
    dist.index = dist.index.astype(int)

    # Fit power law
    x = np.log(dist.index.values).reshape([-1,1])
    y = np.log(dist.values).reshape([-1,1])
    model = LinearRegression()
    model.fit(x,y)

    # Store degree distribution, model, and R2 in uns["scale_free_topology"]
    output = {
        "degree_distribution": dist,
        "model": model,
        "r2": model.score(x,y)
    }

    return output


if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj_path')
    parser.add_argument('--grn_key', default="grn")
    parser.add_argument('--cluster_key', default=None)
    parser.add_argument('--cluster', default=None)

    args = parser.parse_args()

    mdata = mudata.load(args.mudataObj_path)
    _compute_scale_free_topology_fit(mdata, grn_key=args.grn_key, cluster_key=args.cluster_key, cluster=args.cluster)
