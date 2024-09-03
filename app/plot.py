import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import base64
import matplotlib.pyplot as plt
import matplotlib
from io import BytesIO
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
import seaborn as sns



def generate_large_colormap(num_colors):
    """Generate a large colormap using matplotlib

    Parameters
    ----------
    num_colors : int
        Number of colors to generate
    
    Returns
    -------
    list
        List of colors in hex format
    """

    cmap = plt.get_cmap('tab20b', num_colors)
    colors = [matplotlib.colors.rgb2hex(cmap(i)[:3]) for i in range(cmap.N)]
    return colors


def map_categories_to_colors(
    categories
):
    """Map categories to colors
    
    Parameters
    ----------
    categories : list
        List of categories
    
    Returns
    -------
    list
        List of colors in hex format
    dict
        Dictionary mapping categories to colors
    """
    unique_categories = sorted(categories.unique())
    colors = generate_large_colormap(len(unique_categories))
    color_map = {category: color for category, color in zip(unique_categories, colors)}
    return categories.map(color_map).tolist(), color_map


def fig_to_uri(in_fig, close_all=True, **save_args):
    """
    Save a figure as a URI
    
    Parameters
    ----------
    in_fig : matplotlib.figure.Figure
        The input figure
    close_all : bool
        Whether to close all figures
    save_args : dict
        Additional arguments to pass to savefig
    
    Returns
    -------
    str
        URI of the figure
    """
    out_img = BytesIO()
    in_fig.savefig(out_img, format='png', **save_args)
    if close_all:
        in_fig.clf()
        plt.close('all')
    out_img.seek(0)  # rewind file
    encoded = base64.b64encode(out_img.read()).decode("ascii").replace("\n", "")
    return "data:image/png;base64,{}".format(encoded)


def scatterplot(
    data: pd.DataFrame,
    x_column: str,
    y_column: str,
    title: str,
    sorted: bool = True,
    x_axis_title: str = None,
    y_axis_title: str = None,
    cumulative: bool = False,
    show_xaxis_labels: bool = False,
    colors: list = None,  # New parameter for optional colors
    size: int = 1
):
    """Create a scatter plot layout in Dash using Plotly.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame containing the data for the plot.
    x_column : str
        The column to use for the x-axis.
    y_column : str
        The column to use for the y-axis.
    title : str
        Title of the plot.
    x_axis_title : str, optional
        Title for the x-axis.
    y_axis_title : str, optional
        Title for the y-axis.
    cumulative : bool, optional
        Whether to plot cumulative values.
    show_xaxis_labels : bool, optional
        Whether to show labels on the x-axis.
    colors : list, optional
        List of colors corresponding to each point in the plot.

    Returns
    -------
    go.Figure
        A Plotly Figure containing the scatter plot.
    """
    # Compute cumulative values
    if cumulative:
        data[y_column] = data[y_column].cumsum()

    # Sort
    x_data = data.sort_values(y_column, ascending=cumulative)[x_column] if sorted else data[x_column]
    y_data = data.sort_values(y_column, ascending=cumulative)[y_column] if sorted else data[y_column]

    # Plot
    fig = go.Figure(
        data=go.Scattergl(
            x=x_data,
            y=y_data,
            mode='markers',
            marker=dict(
                color=colors if colors is not None else 'blue',  # Use passed colors or default to blue
                size=size
            )
        )
    )

    # Update layout
    fig.update_layout(
        title=title,
        xaxis_title=x_axis_title if x_axis_title else x_column,
        yaxis_title=y_axis_title if y_axis_title else y_column,
        xaxis=dict(showticklabels=show_xaxis_labels),
        plot_bgcolor='rgba(0,0,0,0)',  # Transparent background
    )
    
    return fig


def scatterplot_static(
    ax: plt.Axes,
    data: pd.DataFrame,
    x_column: str,
    y_column: str,
    title: str,
    sorted: bool = True,
    x_axis_title: str = None,
    y_axis_title: str = None,
    cumulative: bool = False,
    show_xaxis_labels: bool = False,
    show_yaxis_labels: bool = False,
    colors: list = None,  # New parameter for optional colors
    cmap: str = None,
    size: int = 1,
):
    """Generate a static scatter plot using Matplotlib."""

    # Give white background
    ax.set_facecolor('white')

    # Compute cumulative values
    if cumulative:
        data[y_column] = data[y_column].cumsum()

    # Sort
    x_data = data.sort_values(y_column, ascending=cumulative)[x_column] if sorted else data[x_column]
    y_data = data.sort_values(y_column, ascending=cumulative)[y_column] if sorted else data[y_column]

    # Plot
    ax.scatter(
        x_data,
        y_data,
        c=colors,
        s=size,
        cmap=cmap
    )

    # remove the top and right spines from plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel(x_axis_title if x_axis_title else x_column)
    ax.set_ylabel(y_axis_title if y_axis_title else y_column)
    ax.set_title(title)
    if not show_xaxis_labels:
        ax.set_xticks([])
    if not show_yaxis_labels:
        ax.set_yticks([])


def barplot(
    data: pd.DataFrame,
    x_column: str,
    y_column: str,
    title: str,
    x_axis_title: str = None,
    y_axis_title: str = None,
    show_xaxis_labels: bool = True
):
    """Create a layout for a filtered bar plot.

    Parameters
    ----------
    app : Dash
        The Dash app instance.
    data : pd.DataFrame
        DataFrame containing the data for the plot.
    x_column : str
        The column to use for the x-axis.
    y_column : str
        The column to use for the y-axis.
    filter_column : str
        The column to apply the filter on.
    starting_filter_value : float
        Initial filter value.
    title : str
        Title of the plot.
    x_axis_title : str, optional
        Title for the x-axis.
    y_axis_title : str, optional
        Title for the y-axis.
    id_suffix : str, optional
        Suffix to add to the id of the input and output components to ensure uniqueness.
    show_xaxis_labels : bool, optional
        Whether to show x-axis labels.

    Returns
    -------
    html.Div
        A Div containing the filter controls and bar plot.
    """
    fig = px.bar(
        data,
        x=x_column,
        y=y_column,
        template='plotly_white'
    ).update_layout(
        title=title,
        xaxis_title=x_axis_title if x_axis_title else x_column,
        yaxis_title=y_axis_title if y_axis_title else y_column,
        xaxis=dict(showticklabels=show_xaxis_labels)
    )
    return fig


def lollipop_plot(
    data: pd.DataFrame,
    x_column: str,
    y_column: str,
    title: str,
    x_axis_title: str = None,
    y_axis_title: str = None,
    show_xaxis_labels: bool = True,
    marker_colors: list = None,
    line_colors: list = None
):
    """Create a lollipop plot with colored markers and lines based on the sign of the y-axis values."""

    fig = go.Figure()

    # Add the vertical lines for the lollipops
    for i, row in data.iterrows():
        fig.add_shape(
            type="line",
            x0=row[x_column],
            y0=0,
            x1=row[x_column],
            y1=row[y_column],
            line=dict(color=line_colors[i])
        )

    # Add the markers for the lollipops
    fig.add_trace(go.Scattergl(
        x=data[x_column],
        y=data[y_column],
        mode='markers',
        marker=dict(color=marker_colors, size=8),
        text=data[x_column]
    ))

    # Update layout
    fig.update_layout(
        title=title,
        xaxis_title=x_axis_title if x_axis_title else x_column,
        yaxis_title=y_axis_title if y_axis_title else y_column,
        xaxis=dict(showticklabels=show_xaxis_labels),
        plot_bgcolor='rgba(0,0,0,0)',  # Transparent background
    )

    return fig


def stacked_barplot(
    data: pd.DataFrame,
    x_column: str,
    y_column: str,
    stack_column: str,
    title: str,
    x_axis_title: str = None,
    y_axis_title: str = None,
    show_xaxis_labels: bool = True
):
    """Create a layout for a filtered stacked bar plot.

    Parameters
    ----------
    app : Dash
        The Dash app instance.
    data : pd.DataFrame
        DataFrame containing the data for the plot.
    x_column : str
        The column to use for the x-axis.
    y_column : str
        The column to use for the y-axis.
    stack_column : str
        The column to use for stacking.
    filter_column : str
        The column to apply the filter on.
    starting_filter_value : float
        Initial filter value.
    title : str
        Title of the plot.
    x_axis_title : str, optional
        Title for the x-axis.
    y_axis_title : str, optional
        Title for the y-axis.
    id_suffix : str, optional
        Suffix to add to the id of the input and output components to ensure uniqueness.
    show_xaxis_labels : bool, optional
        Whether to show x-axis labels.

    Returns
    -------
    html.Div
        A Div containing the filter controls and stacked bar plot.
    """


    fig = px.bar(
        data,
        x=x_column,
        y=y_column,
        color=stack_column,
        barmode='stack',
        template='plotly_white'
    ).update_layout(
        title=title,
        xaxis_title=x_axis_title if x_axis_title else x_column,
        yaxis_title=y_axis_title if y_axis_title else y_column,
        xaxis=dict(showticklabels=show_xaxis_labels)
    )

    return fig


def boxplot(
    data: pd.DataFrame,
    x_column: str,
    y_column: str,
    title: str,
    x_axis_title: str = None,
    y_axis_title: str = None,
    show_xaxis_labels: bool = True
):
    """Create a layout for a box plot.

    Parameters
    ----------
    app : Dash
        The Dash app instance.
    data : pd.DataFrame
        DataFrame containing the data for the plot.
    x_column : str
        The column to use for the x-axis.
    y_column : str
        The column to use for the y-axis.
    filter_column : str
        The column to apply the filter on.
    starting_filter_value : float
        Initial filter value.
    title : str
        Title of the plot.
    x_axis_title : str, optional
        Title for the x-axis.
    y_axis_title : str, optional
        Title for the y-axis.
    id_suffix : str, optional
        Suffix to add to the id of the input and output components to ensure uniqueness.
    show_xaxis_labels : bool, optional
        Whether to show x-axis labels.

    Returns
    -------
    html.Div
        A Div containing the filter controls and box plot.
    """
    fig = px.box(
        data,
        x=x_column,
        y=y_column,
        template='plotly_white'
    ).update_layout(
        title=title,
        xaxis_title=x_axis_title if x_axis_title else x_column,
        yaxis_title=y_axis_title if y_axis_title else y_column,
        xaxis=dict(showticklabels=show_xaxis_labels)
    )
    return fig


def volcano_plot(
    data: pd.DataFrame,
    effect_size_var: str,
    sig_var: str,
    sig_threshold: float,
    effect_size_threshold: float,
    hover_data: list = None
):
    
    # Plot volcano plot
    fig = px.scatter(
        data,
        x=effect_size_var,
        y=sig_var,
        title="",
        hover_data=hover_data
    )

    # Customize layout
    fig.update_layout(
        xaxis_title=effect_size_var,
        yaxis_title=sig_var,
        yaxis=dict(tickformat=".1f"),
        width=1000,
        height=800,
        xaxis_tickfont=dict(size=4),
        plot_bgcolor='rgba(0,0,0,0)',  # Transparent background
    )

    # Add horizontal dashed line for significance threshold (make it red)
    fig.add_hline(
        y=-np.log10(sig_threshold),
        line_dash="dash",
        line_color="red",
        annotation_text=f'Significance Threshold ({sig_threshold})',
        annotation_position="top right"
    )

    # Add vertical dashed lines for effect size threshold
    fig.add_vline(
        x=effect_size_threshold,
        line_dash="dash",
        line_color="red",
        annotation_text=f'Effect Size Threshold ({effect_size_threshold})',
        annotation_position="top right"
    )

    return fig


def convertDatTraits(data):
    """
    get data trait module base on samples information

    :return: a dataframe contains information in suitable format for plotting module trait relationship heatmap
    :rtype: pandas dataframe
    """
    datTraits = pd.DataFrame(index=data.index)
    for i in range(data.shape[1]):
        data.iloc[:, i] = data.iloc[:, i].astype(str)
        if len(np.unique(data.iloc[:, i])) == 2:
            datTraits[data.columns[i]] = data.iloc[:, i]
            org = np.unique(data.iloc[:, i]).tolist()
            rep = list(range(len(org)))
            datTraits.replace(to_replace=org, value=rep,
                              inplace=True)
        elif len(np.unique(data.iloc[:, i])) > 2:
            for name in np.unique(data.iloc[:, i]):
                datTraits[name] = data.iloc[:, i]
                org = np.unique(data.iloc[:, i])
                rep = np.repeat(0, len(org))
                rep[np.where(org == name)] = 1
                org = org.tolist()
                rep = rep.tolist()
                datTraits.replace(to_replace=org, value=rep, inplace=True)

    return datTraits


def plot_topic_trait_relationship_heatmap(
        ax,
        cell_topic_participation,
        metaData,
        covariates,
    ):
    """
    plot topic-trait relationship heatmap
    """
    datTraits = convertDatTraits(metaData[covariates])
    datTraits.index = cell_topic_participation.index

    topicsTraitCor = pd.DataFrame(index=cell_topic_participation.columns,
                                  columns=datTraits.columns,
                                  dtype="float")
    topicsTraitPvalue = pd.DataFrame(index=cell_topic_participation.columns,
                                     columns=datTraits.columns,
                                     dtype="float")
    min_cell_participation = cell_topic_participation.min().min()
    for i in cell_topic_participation.columns:
        for j in datTraits.columns:
            tmp = cell_topic_participation[
                ~np.isclose(cell_topic_participation[i],
                            min_cell_participation, atol=min_cell_participation)]
            tmp = stats.spearmanr(tmp[i], datTraits.loc[tmp.index, j], alternative='greater')
            topicsTraitCor.loc[i, j] = tmp[0]
            topicsTraitPvalue.loc[i, j] = tmp[1]

    topicsTraitCor.fillna(0.0, inplace=True)
    topicsTraitPvalue.fillna(1.0, inplace=True)

    for i in range(topicsTraitPvalue.shape[0]):
        rejected, tmp = fdrcorrection(topicsTraitPvalue.iloc[i, :])
        if not rejected.all():
            topicsTraitPvalue.iloc[i, :] = tmp

    xlabels = cell_topic_participation.columns
    ylabels = datTraits.columns

    sns.set(font_scale=1.5)
    res = sns.heatmap(topicsTraitCor.T, cmap='RdBu_r',
                        vmin=-1, vmax=1, ax=ax, annot_kws={'size': 12, "weight": "bold"},
                        xticklabels=xlabels, yticklabels=ylabels)
    res.set_xticklabels(res.get_xmajorticklabels(), fontsize=12, fontweight="bold", rotation=90)
    res.set_yticklabels(res.get_ymajorticklabels(), fontsize=12, fontweight="bold")
    plt.yticks(rotation=0)
    ax.set_facecolor('white')
