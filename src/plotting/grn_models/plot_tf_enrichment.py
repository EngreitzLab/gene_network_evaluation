import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
import scanpy as sc


def plot_tf_activity_vs_expression(
    tf_acts,
    tf_expr,
    tf,
    layer=None,
):
    with sns.plotting_context("notebook"):
        fig, ax = plt.subplots(1, 2, figsize=(8, 4))
        sc.pl.violin(tf_acts, groupby="annotation", keys=tf, show=False, ax=ax[0])
        ax[0].set_xticklabels(ax[0].get_xticklabels(), rotation=90, ha="right")
        ax[0].set_ylabel(f"{tf} activity")
        sc.pl.violin(tf_expr, groupby="annotation", keys=tf, show=False, ax=ax[1], layer=layer, use_raw=False)
        ax[1].set_xticklabels(ax[1].get_xticklabels(), rotation=90, ha="right")
        ax[1].set_ylabel(f"{tf} expression")
        plt.tight_layout()

def plot_tf_activity_on_dim_reduction(
    tf_acts,
    dim_reduction_key,
    tf,
    cluster_key=None,
):
    with sns.plotting_context("notebook"):
        if cluster_key is None:
            _, ax = plt.subplots(1, 2, figsize=(8, 4))
        else:
            _, ax = plt.subplots(1, 3, figsize=(12, 3.5))
        sc.pl.violin(tf_acts, groupby=cluster_key, keys=tf, show=False, ax=ax[0])
        sc.pl.embedding(tf_acts, basis=dim_reduction_key, color=[tf], s=25, vmax="p99.9", cmap="viridis", ax=ax[1], show=False)
        if cluster_key is not None:
            sc.pl.embedding(tf_acts, basis=dim_reduction_key, color=[cluster_key], s=25, ax=ax[2], show=False)
        for tick in ax[0].get_xticklabels():
            tick.set_rotation(90)
        plt.tight_layout()
        ax[0].set_title("")