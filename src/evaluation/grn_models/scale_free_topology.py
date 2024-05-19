import os
import argparse

import mudata

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from igraph import Graph
from sklearn.linear_model import LinearRegression


def compute_scale_free_topology_fit(
    mdata: mudata.MuData,
    grn_key: str = "grn",
    cluster_key: str = None,
    cluster: str = None,
    inplace: bool = True
):
    # Extract GRN
    grn = mdata.uns[grn_key].copy()
    if cluster_key is not None:
        assert cluster is not None, "cluster must be provided if cluster_key is provided"
        grn = grn.query("cluster_key == @cluster")

    # Build iGraph
    g = Graph.DataFrame(grn[["source", "target"]], directed=True, use_vids=False)
    g.es["weight"] = grn["weight"].values.copy()

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

    # Update mdata
    if inplace:
        mdata.uns["scale_free_topology"] = output
    else:
        return output


def plot_scale_free_topology_fit(
    mdata: mudata.MuData,
    uns_key: str = "scale_free_topology",
    ax=None
):
    # Create figure
    fig, ax = plt.subplots(1, 2, figsize=(8,4))

    # Extract degree distribution and model
    dist = mdata.uns[uns_key]["degree_distribution"]
    model = mdata.uns[uns_key]["model"]

    # Plot degree (k) vs P(k)
    ax[0].scatter(dist.index.values, dist.values, c="black")
    ax[0].set_title("degree distribution")
    ax[0].set_xlabel("k")
    ax[0].set_ylabel("P(k)")

    # Get log-log scale
    x = np.log(dist.index.values).reshape([-1,1])
    y = np.log(dist.values).reshape([-1,1])

    # For each x value, calculate the predicted y value
    x_ = np.linspace(x.min(), x.max(), 100).reshape([-1,1])
    y_ = model.predict(x_)

    # Plot log-log scale
    ax[1].set_title(f"degree distribution (log scale)\nslope: {model.coef_[0][0] :.4g}, r2: {model.score(x,y) :.4g}")
    ax[1].plot(x_.flatten(), y_.flatten(), c="black", alpha=0.5)
    ax[1].scatter(x.flatten(), y.flatten(), c="black")
    ax[1].set_ylim([y.min()-0.2, y.max()+0.2])
    ax[1].set_xlim([-0.2, x.max()+0.2])
    ax[1].set_xlabel("log k")
    ax[1].set_ylabel("log P(k)")

    # Return figure
    return fig, ax


if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj_path')
    parser.add_argument('--grn_key', default="grn")
    parser.add_argument('--cluster_key', default=None)
    parser.add_argument('--cluster', default=None)

    args = parser.parse_args()

    mdata = mudata.load(args.mudataObj_path)
    compute_scale_free_topology_fit(mdata, grn_key=args.grn_key, cluster_key=args.cluster_key, cluster=args.cluster)
