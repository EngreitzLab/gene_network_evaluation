import pandas as pd
from dash import html, dcc
from dash.dash_table import DataTable
from dash import Input, Output
import plotly.express as px
import plotly.io as pio
from data_processing import count_unique, filter_data

# Set the default template to plotly_white
pio.templates.default = "plotly_white"

def create_scatter_layout(
    data: pd.DataFrame,
    x_column: str,
    y_column: str,
    title: str,
    sorted: bool = True,
    x_axis_title: str = None,
    y_axis_title: str = None,
    cumulative: bool = False,
    id_suffix: str = '',
    show_xaxis_labels: bool = False
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

    Returns
    -------
    html.Div
        A Div containing the scatter plot.
    """
    # Compute cumulative values
    if cumulative:
        data[y_column] = data[y_column].cumsum()

    # Sort
    x_data = data.sort_values(y_column, ascending=cumulative)[x_column] if sorted else data[x_column]
    y_data = data.sort_values(y_column, ascending=cumulative)[y_column] if sorted else data[y_column]
    
    # Plot
    fig = px.scatter(
            x=x_data,
            y=y_data,
            template='plotly_white'
        ).update_layout(
            title=title, 
            xaxis_title=x_axis_title if x_axis_title else x_column, 
            yaxis_title=y_axis_title if y_axis_title else y_column,
            xaxis=dict(showticklabels=show_xaxis_labels)
        )
    

    # Update hover template to show scientific notation
    fig.update_traces(hovertemplate='%{y:e}')

    # Return the layout
    return html.Div(
        [dcc.Graph(figure=fig)],
        id=f'{x_column}-scatter-container-{id_suffix}'
    )


def create_filtered_barplot_layout(
    app,
    data: pd.DataFrame,
    x_column: str,
    y_column: str,
    filter_column: str,
    starting_filter_value: float,
    title: str,
    x_axis_title: str = None,
    y_axis_title: str = None,
    id_suffix: str = '',
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
    input_id = f'{filter_column}-input-{id_suffix}'
    container_id = f'{filter_column}-barplot-container-{id_suffix}'
    print(f"Creating layout for input id: {input_id} and container id: {container_id}")
    layout = html.Div([
        html.Label(f'Filter by {filter_column}'),
        dcc.Input(
            id=input_id,
            type='number',
            value=starting_filter_value,
        ),
        html.Div(id=container_id)
    ])

    @app.callback(
        Output(container_id, 'children'),
        Input(input_id, 'value')
    )
    def update_barplot(filter_value):
        print("Callback triggered for filter value:", filter_value)
        filtered_data = data[data[filter_column] <= filter_value]
        if filtered_data.empty:
            return html.Div("No data available for the selected filter.")
        unique_df = count_unique(categorical_var=x_column, count_var=y_column, dataframe=filtered_data)
        fig = px.bar(
            unique_df.sort_values(y_column, ascending=False),
            x=x_column,
            y=y_column,
            template='plotly_white'
        ).update_layout(
            title=title,
            xaxis_title=x_axis_title if x_axis_title else x_column,
            yaxis_title=y_axis_title if y_axis_title else y_column,
            xaxis=dict(showticklabels=show_xaxis_labels)
        )
        return dcc.Graph(figure=fig)

    return layout


def create_filtered_stacked_barplot_layout(
    app,
    data: pd.DataFrame,
    x_column: str,
    y_column: str,
    stack_column: str,
    filter_column: str,
    starting_filter_value: float,
    title: str,
    x_axis_title: str = None,
    y_axis_title: str = None,
    id_suffix: str = '',
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
    input_id = f'{filter_column}-input-{id_suffix}'
    container_id = f'{filter_column}-stacked-barplot-container-{id_suffix}'
    print(f"Creating layout for input id: {input_id} and container id: {container_id}")

    layout = html.Div([
        html.Label(f'Filter by {filter_column}'),
        dcc.Input(
            id=input_id,
            type='number',
            value=starting_filter_value,
        ),
        html.Div(id=container_id)
    ])

    @app.callback(
        Output(container_id, 'children'),
        Input(input_id, 'value')
    )
    def update_stacked_barplot(filter_value):
        print("Callback triggered for filter value:", filter_value)
        # Filter data
        filtered_data = data[data[filter_column] <= filter_value]

        # Compute counts for each stack value
        stacked_data_list = []
        for stack_value in filtered_data[stack_column].unique():
            subset = filtered_data[filtered_data[stack_column] == stack_value]
            unique_counts = count_unique(categorical_var=x_column, count_var=y_column, dataframe=subset)
            unique_counts[stack_column] = stack_value
            stacked_data_list.append(unique_counts)

        # Handle empty data case
        if len(stacked_data_list) == 0:
            return html.Div("No data available for the selected filter.")

        # Combine stacked data
        stacked_data = pd.concat(stacked_data_list)

        # Sort by total counts
        total_counts = stacked_data.groupby(x_column)[y_column].sum().reset_index()
        total_counts = total_counts.rename(columns={y_column: "Total"})
        stacked_data = stacked_data.merge(total_counts, on=x_column)
        stacked_data = stacked_data.sort_values("Total", ascending=False)

        fig = px.bar(
            stacked_data,
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

        return dcc.Graph(figure=fig)

    return layout
