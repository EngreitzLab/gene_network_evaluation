import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import io
import base64
import matplotlib.pyplot as plt
import matplotlib


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
    size: int = 8
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
    save_as_base64: bool = False
):
    """Generate a static scatter plot using Matplotlib."""

    # Compute cumulative values
    if cumulative:
        data[y_column] = data[y_column].cumsum()

    # Sort
    x_data = data.sort_values(y_column, ascending=cumulative)[x_column] if sorted else data[x_column]
    y_data = data.sort_values(y_column, ascending=cumulative)[y_column] if sorted else data[y_column]

    # Plot
    fig, ax = plt.subplots()
    scatter = ax.scatter(
        x_data,
        y_data,
        c=colors,
        s=size,
        cmap=cmap
    )

    if cmap:
        fig.colorbar(scatter, ax=ax)

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
    
    # Save the plot to a BytesIO object
    if save_as_base64:
        buf = io.BytesIO()
        plt.savefig(buf, format='png')
        plt.close(fig)
        buf.seek(0)
        
        # Encode the plot as a base64 string
        img_base64 = base64.b64encode(buf.read()).decode('ascii')
        
        return f"data:image/png;base64,{img_base64}"
    
    else:
        return fig
    

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
