import os
import argparse

import mudata

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def plot_scale_free_topology_fit(
    degree_dist,
    model,
    ax=None
):
    # Create figure
    fig, ax = plt.subplots(1, 2, figsize=(8,4))

    # Plot degree (k) vs P(k)
    ax[0].scatter(degree_dist.index.values, degree_dist.values, c="black")
    ax[0].set_title("degree distribution")
    ax[0].set_xlabel("k")
    ax[0].set_ylabel("P(k)")

    # Get log-log scale
    x = np.log(degree_dist.index.values).reshape([-1,1])
    y = np.log(degree_dist.values).reshape([-1,1])

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


def plot_goodness_of_fit(
    fit_df
):
    sns.jointplot(
        data=fit_df,
        x="r2",
        y="n_regulators",
        kind="scatter",
        marginal_kws=dict(bins=30, fill=True),
        marginal_ticks=True,
        s=10,
        alpha=0.5,
    )