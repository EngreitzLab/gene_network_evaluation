import dash
import dash_bootstrap_components as dbc
from dash import html, dcc, callback, Input, Output, dash_table
import plotly.graph_objects as go
import diskcache
from plot import map_categories_to_colors, scatterplot, scatterplot_static
import matplotlib.pyplot as plt
import yaml

# Register the page
dash.register_page(__name__, order=0, path='/')

# Get the app and cache
app = dash.get_app()
cache = diskcache.Cache("./.cache")
results = cache.get("results")

# Grab the dimensionality reduction data and categorical keys
obsms = results["obsms"]
default_obsm = 'X_umap' if 'X_umap' in obsms else list(obsms.keys())[0]

# Get the columns that don't have the selected dim reduction prefix
categorical_keys = results["categorical_keys"]
continuous_keys = results["continuous_keys"]
covariate_keys = categorical_keys + continuous_keys
default_covariate = categorical_keys[0] if categorical_keys else None

# Generate colors for the categorical keys
categorical_key_colors = {}
categorical_key_color_maps = {}
for key in categorical_keys:
    colors, color_map = map_categories_to_colors(results["obs"][key])
    categorical_key_colors[key] = colors
    categorical_key_color_maps[key] = color_map

# Get the path evaluation config
evaluation_config = results["evaluation_config"]
evaluation_yaml = yaml.dump(evaluation_config, default_flow_style=False, sort_keys=False)

# Get the software versions
software_versions = results["software_versions"]
software_versions_yaml = yaml.dump(software_versions, default_flow_style=False, sort_keys=False)

# Create the layout
layout = dbc.Container([
    html.H1("Data Overview", className="mb-4"),

    # Covariate and obsm dropdowns at the top
    dbc.Row([

        # Covariate dropdown
        dbc.Col([
            html.Label("Select Covariate"),
            dcc.Dropdown(
                id='covariate-selector',
                options=[{'label': key, 'value': key} for key in covariate_keys],
                value=default_covariate,
                clearable=False
            )
        ], width=6),

        # Obsm dropdown
        dbc.Col([
            html.Label("Select Dimensionality Reduction"),
            dcc.Dropdown(
                id='dim-reduction-selector',
                options=[{'label': key, 'value': key} for key in obsms.keys()],
                value=default_obsm,
                clearable=False
            )
        ], width=6),
    ], className="mb-4"),

    # Main
    dbc.Row([
        
        # Dimensionality Reduction Plot
        dbc.Col([
            dcc.Graph(id='dim-reduction-plot', className="mt-3"),
        ], width=6),  # Larger area for plot
        
        # Legend occupies some
        dbc.Col([
            html.Div(id='legend', className="mt-3")
        ], width=2),
        
        # Barplot/Histogram
        dbc.Col([
            dcc.Graph(id='covariate-distribution', className="mt-3"),
        ], width=4),

    ], className="mb-4"),

    # Summary of evaluation parameters
    dbc.Row([
        dbc.Col([
            html.H3("Evaluation Parameters", className="mt-4"),
            dcc.Markdown(
                f"```yaml\n{evaluation_yaml}\n```", 
                style={
                    'whiteSpace': 'pre-wrap', 
                    'backgroundColor': '#f8f9fa', 
                    'padding': '20px',
                    'borderRadius': '10px',
                    'border': '1px solid #dee2e6',
                    'fontFamily': 'Courier, monospace'
                }
            )
        ], width=12),
    ], className="mb-4"),

    # Software versions
    dbc.Row([
        dbc.Col([
            html.H3("Package versions for evaluation", className="mt-4"),
            dcc.Markdown(
                f"```yaml\n{software_versions_yaml}\n```",
                style={
                    'whiteSpace': 'pre-wrap',
                    'backgroundColor': '#f8f9fa',
                    'padding': '20px',
                    'borderRadius': '10px',
                    'border': '1px solid #dee2e6',
                    'fontFamily': 'Courier, monospace'
                }
            )
        ], width=12),
    ]),
], fluid=True, className="p-4")


# Callback to update the scatter plot and legend based on selected dim reduction and covariate
@callback(
    [Output('dim-reduction-plot', 'figure'),
     Output('legend', 'children')],
    [Input('dim-reduction-selector', 'value'), 
     Input('covariate-selector', 'value')]
)
def update_dim_reduction_plot(
    selected_dim_reduction, 
    selected_covariate,
    size=1,
    debug=True,
    static=True
):

    if debug:
        print(f"Selected dim reduction: {selected_dim_reduction}, selected covariate: {selected_covariate}")
        print(f"Static: {static}")

    # Check if the selected dim reduction and covariate are valid
    if selected_dim_reduction is None or selected_covariate is None:
        return go.Figure(), html.Div()  # Return an empty figure and legend placeholder

    # grab cell metadata
    cell_metadata = results["obs"]

    # Extract the relevant data for plotting (they will be selected_dim_reduction_0 and selected_dim_reduction_1)
    data = obsms[selected_dim_reduction][[f'{selected_dim_reduction}_0', f'{selected_dim_reduction}_1']].loc[cell_metadata["barcode"]]
    data.columns = ['X', 'Y']

    # Add the selected covariate to the data
    data[selected_covariate] = cell_metadata[selected_covariate].values

    # Get covariate categories to more colors
    if selected_covariate in categorical_keys:
        colors = categorical_key_colors[selected_covariate]
        categorical_color_map = categorical_key_color_maps[selected_covariate]
        continuous_color_map = None
        if debug:
            print(f"Categorical covariate: {selected_covariate}")
            print(f"Colors: {colors[:5]}")
            print(f"Color map: {categorical_color_map}")
    else:
        colors = data[selected_covariate]
        categorical_color_map = None
        continuous_color_map = None
        if debug:
            print(f"Continuous covariate: {selected_covariate}")
            print(f"Colors: {colors[:5]}")
            print(f"Color map: {continuous_color_map}")

    if static:

        # Plot using static matplotlib function
        fig = scatterplot_static(
            data=data,
            x_column="X",
            y_column="Y",
            sorted=False,
            title='',
            x_axis_title=f'{selected_dim_reduction}1',
            y_axis_title=f'{selected_dim_reduction}2',
            colors=colors,
            cmap=continuous_color_map,
            size=size,
        )
        
        # add legend to the plot
        if selected_covariate in categorical_keys:
            fig.legend(handles=[plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=label) for label, color in categorical_color_map.items()])

        # Return the figure and an empty Div
        return fig, html.Div()
    
    else:
        
        # Plot using dynamic plotly function
        fig = scatterplot(
            data=data,
            x_column="X",
            y_column="Y",
            sorted=False,
            title='',
            x_axis_title=f'{selected_dim_reduction}1',
            y_axis_title=f'{selected_dim_reduction}2',
            colors=colors,
            size=1
        )

        # Generate the legend
        if selected_covariate in categorical_keys:
            legend_items = [html.Div([
                                html.Span(style={'backgroundColor': color, 'display': 'inline-block', 'width': '20px', 'height': '20px', 'marginRight': '10px'}),
                                html.Span(category)
                            ], style={'marginBottom': '5px'}) for category, color in color_map.items()]
            legend = html.Div(legend_items, style={'display': 'flex', 'flexDirection': 'column'})
        else:
            legend = html.Div()
        
        # Return the figure and legend
        return fig, legend


# Callback to update the barplot based on selected covariate
@callback(
    Output('covariate-distribution', 'figure'),
    [Input('covariate-selector', 'value')]
)
def update_covariate_distribution_plot(
    selected_covariate,
    debug=False
):

    if debug:
        print(f"Selected covariate: {selected_covariate}")

    # Check if the selected covariate is valid
    if selected_covariate is None:
        return go.Figure()  # Return an empty figure
    
    # grab cell metadata
    cell_metadata = results["obs"]
    
    # Extract the relevant data for plotting
    data = cell_metadata[selected_covariate]
    
    # if the selected covariate is continuous, plot a histogram
    if selected_covariate in continuous_keys:
        # Plot a histogram of the selected covariate with white background
        fig = go.Figure(
            go.Histogram(
                x=data,
                histnorm='density',
                marker_color='blue'
            )
        )
        fig.update_layout(
            title=f"Distribution of {selected_covariate}",
            xaxis_title=selected_covariate,
            yaxis_title="Density",
            plot_bgcolor='white',
            showlegend=False
        )
    
    # if the selected covariate is categorical, plot a count barplot with white background
    else:
        color_map = categorical_key_color_maps[selected_covariate]
        counts = data.value_counts()
        categories = counts.index
        colors = [color_map[category] for category in categories]
        fig = go.Figure(
            go.Bar(
                x=categories,
                y=counts,
                marker_color=colors
            )
        )

        fig.update_layout(
            title=f"Distribution of {selected_covariate}",
            xaxis_title=selected_covariate,
            yaxis_title="Count",
            plot_bgcolor='white',
            showlegend=False
        )

    return fig
