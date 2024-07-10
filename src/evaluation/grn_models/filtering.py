import numpy as np
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO)


def filter_network(
    network: pd.DataFrame,
    source: str = "tf",
    target: str = "gene",
    pval: float = 0.001,
    score: float = 0,
    ntop: int = 2000,
    cluster_key: str = None,
    cluster: str = None,
    min_regulators: int = 0,
    min_targets: int = 0,
):
    """
    Filter a gene regulatory network (filtered_network) based on specified criteria.

    Parameters
    ----------
    filtered_network : pandas.DataFrame
        The filtered_network as a DataFrame.
    pval : float, optional
        P-value threshold for filtering edges, by default 0.001.
    score : float, optional
        Score threshold for filtering edges, by default 0.
    ntop : int, optional
        Number of top regulators to keep, by default 2000.
    cluster_key : str, optional
        Key for the cluster column in the filtered_network DataFrame, by default None.
    cluster : str, optional
        The cluster of interest to filter the filtered_network for, by default None.
    min_regulators : int, optional
        Minimum number of regulators a gene must have, by default 0.
    min_targets : int, optional
        Minimum number of targets a regulator must have, by default 0.

    Returns
    -------
    pandas.DataFrame
        The filtered filtered_network dataframe.

    """
    logger = logging.getLogger(__name__)

    # Log initial number of edges
    logger.info(f"Initial network has {len(network)} edges")
    filtered_network = network.copy()
    
    # Choose a cluster of interest
    if cluster_key is not None:
        if cluster is None or cluster_key not in filtered_network.columns:
            raise ValueError("If cluster_key is provided, cluster must also be provided and cluster_key must be in filtered_network.columns")
        filtered_network = filtered_network[filtered_network[cluster_key] == cluster]
    logger.info(f"filtered_network for {cluster} has {len(filtered_network)} edges")
    
    # Filter based on pval and score thresholds
    if pval is not None:
        filtered_network = filtered_network[filtered_network["pval"] <= pval]
    if score is not None:
        filtered_network = filtered_network[np.abs(filtered_network["score"]) >= score]
    logger.info(f"filtered_network after filtering on edge strength has {len(filtered_network)} edges")
    
    # Keep only the top ntop regulators
    if ntop is not None:
        filtered_network["score_abs"] = np.abs(filtered_network["score"])
        filtered_network = filtered_network.sort_values("score_abs", ascending=False)
        filtered_network = filtered_network.head(ntop)
        filtered_network = filtered_network.drop(columns="score_abs")
    logger.info(f"filtered_network after filtering on top regulators has {len(filtered_network)} edges")
    
    # Filter out sources that regulate fewer than min_targets targets
    sources = filtered_network[source].value_counts()
    sources = sources[sources > min_targets].index
    filtered_network = filtered_network[filtered_network[source].isin(sources)]
    logger.info(f"filtered_network after filtering on minimum number of targets has {len(filtered_network)} edges")
    
    # Filter out targets regulated by fewer than min_regulators sources
    targets = filtered_network[target].value_counts()
    targets = targets[targets >= min_regulators].index
    filtered_network = filtered_network[filtered_network[target].isin(targets)]
    logger.info(f"filtered_network after filtering on minimum number of regulators has {len(filtered_network)} edges")
    
    return filtered_network
