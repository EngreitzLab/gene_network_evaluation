import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import seaborn as sns


def plot_module(
    network_data: pd.DataFrame, 
    source_name: str = None,
    target_name: str = None,
    source_col: str = 'source', 
    target_col: str = 'target', 
    score_col: str = 'score', 
    n_top: int = None,
    source_color: str = 'skyblue', 
    target_color: str = 'lightgrey',
    source_shape: str = 's',  # Square by default
    target_shape: str = 'o',  # Circle by default
    cmap=None
) -> None:
    """
    Visualizes a directed network module for a specified source or target node.
    
    Parameters:
    - network_data: DataFrame containing the network data.
    - source_name: The name of the source node (for source-to-target visualization).
    - target_name: The name of the target node (for target-to-source visualization).
    - source_col: Column name for source nodes.
    - target_col: Column name for target nodes.
    - score_col: Column name for interaction scores.
    - source_color: Color for the source node.
    - target_color: Color for the target nodes.
    - source_shape: Shape for the source node.
    - target_shape: Shape for the target nodes.
    - cmap: Color map for edge coloring based on weights.
    """
    # Determine the direction of visualization
    if source_name and not target_name:
        filtered_data = network_data[network_data[source_col] == source_name]
        if n_top is not None:
            filtered_data = filtered_data.nlargest(n_top, score_col)
        main_node = source_name
        secondary_nodes = filtered_data[target_col].unique()
        main_shape = source_shape
        main_color = source_color
        secondary_shape = target_shape
        secondary_color = target_color
        layer_positions = {main_node: (0.5, 1), **{node: (i / max(1, len(secondary_nodes) - 1), 0) for i, node in enumerate(secondary_nodes)}}
    elif target_name and not source_name:
        filtered_data = network_data[network_data[target_col] == target_name]
        if n_top is not None:
            filtered_data = filtered_data.nlargest(n_top, score_col)
        main_node = target_name
        secondary_nodes = filtered_data[source_col].unique()
        main_shape = target_shape
        main_color = target_color
        secondary_shape = source_shape
        secondary_color = source_color
        layer_positions = {main_node: (0.5, 0), **{node: (i / max(1, len(secondary_nodes) - 1), 1) for i, node in enumerate(secondary_nodes)}}
    else:
        raise ValueError("Please specify either a source_name or a target_name, but not both.")

    # Create the directed graph
    G = nx.DiGraph()
    G.add_node(main_node, shape=main_shape, color=main_color, label=main_node)
    weights = []
    for node in secondary_nodes:
        weight = filtered_data[filtered_data[target_col] == node][score_col].values[0] if source_name else filtered_data[filtered_data[source_col] == node][score_col].values[0]
        weights.append(weight)
        G.add_node(node, shape=secondary_shape, color=secondary_color, label=node)
        G.add_edge(main_node if source_name else node, node if source_name else main_node, weight=weight)

    # Determine the appropriate color map if not provided
    if cmap is None:
        if any(weight < 0 for weight in weights) and any(weight > 0 for weight in weights):
            cmap = plt.cm.bwr  # Use blue-white-red for both positive and negative values
        else:
            cmap = plt.cm.Blues if all(weight >= 0 for weight in weights) else plt.cm.Reds  # Use Blues or Reds for positive or negative only

    # Normalize the color scale
    norm = plt.Normalize(min(weights), max(weights))
    edge_colors = [cmap(norm(weight)) for weight in weights]

    # Draw the graph
    plt.figure(figsize=(10, 5))
    pos = layer_positions
    # Draw nodes by shape
    for shape in set([main_shape, secondary_shape]):
        specific_nodes = [node for node, data in G.nodes(data=True) if data['shape'] == shape]
        nx.draw_networkx_nodes(G, pos, nodelist=specific_nodes, node_shape=shape, node_color=[G.nodes[node]['color'] for node in specific_nodes],
                               edgecolors='black', node_size=700)
    nx.draw_networkx_labels(G, pos)
    edges = G.edges(data=True)
    nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color=edge_colors, width=2, arrowstyle='-|>', arrowsize=20)

    # Add a color bar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, label='Interaction Strength', format='%.3f')
    plt.axis('off')  # Turn off the axis
    plt.show()


def plot_grn_statistics(
    network, 
    source_col, 
    target_col, 
    score_col
):
    """
    Generate plots to visualize the GRN statistics.

    Parameters
    ----------
    network : pandas.DataFrame
        DataFrame containing the GRN data.
    source_col : str
        Column name for the source (e.g., transcription factor).
    target_col : str
        Column name for the target (e.g., gene or CRE).
    score_col : str
        Column name for the score (e.g., binding score).

    Returns
    -------
    None
        Displays the plots.

    """
    # Calculate the number of connections per source and target
    num_connections_per_source = network.groupby(source_col).size()
    num_connections_per_target = network.groupby(target_col).size()
    
    with sns.plotting_context("notebook", font_scale=1.5):
        _, ax = plt.subplots(1, 3, figsize=(12, 4))

        # Plot score distribution
        sns.histplot(network[score_col], bins=20, ax=ax[0], log_scale=True)
        ax[0].set_xlabel("Score")
        ax[0].set_ylabel("Num links")
        ax[0].set_title(f"{score_col} distribution")

        sns.violinplot(num_connections_per_target, inner=None, ax=ax[1])
        sns.boxplot(num_connections_per_target, color="white", width=0.2, ax=ax[1])
        ax[1].set_ylabel(f"Num {source_col}")
        ax[1].set_title(f"Num {source_col}s per {target_col}")

        sns.violinplot(num_connections_per_source, inner=None, ax=ax[2])
        sns.boxplot(num_connections_per_source, color="white", width=0.2, ax=ax[2])
        ax[2].set_ylabel(f"Num {target_col}s")
        ax[2].set_title(f"Num {target_col}s per {source_col}")

        plt.tight_layout()
