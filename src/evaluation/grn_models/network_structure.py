import numpy as np
from igraph import Graph

def create_igraph_from_df(
    df,
    source_col="source",
    target_col="target",
    score_col="score",
    directed=True
):
    g = Graph.DataFrame(df[[source_col, target_col]], directed=directed, use_vids=False)
    g.es["weight"] = np.abs(df[score_col].values)
    return g

def node_scores(g):
    df = g.get_vertex_dataframe()
    for i in ["all", "in", "out"]:
        df[f"degree_{i}"] = g.degree(mode=i)
        df[f"degree_centrality_{i}"] = df[f"degree_{i}"] / (df.shape[0]-1)
    df["betweenness_centrality"] = g.betweenness(directed=True, weights="weight")
    df["eigenvector_centrality"] = g.eigenvector_centrality(directed=False, weights="weight")
    df["page_rank"] = g.pagerank(weights="weight")
    df = df.set_index("name")
    df.index.name = None
    return df