import dash
import dash_dangerously_set_inner_html
import dash_bootstrap_components as dbc
from dash import html, dcc, callback, Input, Output
import plotly.graph_objects as go
import mudata
import os
import pandas as pd
from plot import scatterplot

# Ouput directory
path_pipeline_outs = "/cellar/users/aklie/opt/gene_program_evaluation/dashapp/example_data/iPSC_EC_evaluations"
data_key = "rna"

dash.register_page(__name__, order=0)

# Grab html rep from dcc.Store
mudata.set_options(display_style="html", display_html_expand=0b000)
html_rep = dcc.Store(id="mdata", data="")

# Load in mudata html rep from ./mudata.html
with open(os.path.join(path_pipeline_outs, "mudata.html"), "r") as f:
    html_rep = f.read()

# Load in dim reduction .tsvs if it exists (for all files with .tsv extension)
dim_reduce = {}
for file in os.listdir(path_pipeline_outs):
    if file.endswith(".tsv"):
        obsm_key = file.split(".")[0]
        
        # Load in tsv
        df = pd.read_csv(os.path.join(path_pipeline_outs, file), sep="\t")
        
        # Store in dim_reduce
        dim_reduce[obsm_key] = df

# Get the columns that don't have the selected dim reduction prefix
default_dim_reduction = 'X_umap' if 'X_umap' in dim_reduce else list(dim_reduce.keys())[0]
no_reduce_columns = [col for col in dim_reduce[default_dim_reduction].columns if not col.startswith(default_dim_reduction)]
default_covariate = 'rna:leiden' if 'rna:leiden' in dim_reduce[default_dim_reduction] else no_reduce_columns[0]

# Create the layout using tabs for each section
layout = dbc.Container([
    html.H1("Overview", className="mb-4"),

    dcc.Tabs([
        # Tab 1: mdata HTML Representation
        dcc.Tab(label='mdata Overview', children=[
            html.H2("mdata Overview", className="mt-4 mb-3"),
            dash_dangerously_set_inner_html.DangerouslySetInnerHTML(html_rep),
            html.Div(id='mdata-html', className='section', style={'width': '75%', 'display': 'inline-block', 'vertical-align': 'top'}),
        ]),

        # Tab 2: Dimensionality Reduction Visualization
        dcc.Tab(label='Dimensionality Reduction', children=[
            html.H2("Dimensionality Reduction", className="mt-4 mb-3"),
            # Dropdowns for selecting dim reduction and covariate
            dbc.Row([
                dbc.Col([
                    html.Label("Select Dimensionality Reduction"),
                    dcc.Dropdown(
                        id='dim-reduction-selector',
                        options=[{'label': key, 'value': key} for key in dim_reduce.keys()],
                        value=default_dim_reduction,
                        clearable=False
                    )
                ], width=6),
                dbc.Col([
                    html.Label("Select Covariate"),
                    dcc.Dropdown(
                        id='covariate-selector',
                        options=[],  # Will be populated based on dim reduction
                        value=default_covariate,
                        clearable=False
                    )
                ], width=6),
            ], className="mb-4"),

            # Scatter plot and legend side by side
            dbc.Row([
                # Plot occupies half the space
                dbc.Col([
                    dcc.Graph(id='dim-reduction-scatter', className="mt-3"),
                ], width=6),  # Half the width
                
                # Legend occupies the other half
                dbc.Col([
                    html.Div(id='legend', className="mt-3")
                ], width=6)  # Half the width
            ]),
        ]),

        # Tab 3: Run History and Package Versions
        dcc.Tab(label='Run History', children=[
            html.H2("Run History and Package Versions", className="mt-4 mb-3"),
            html.Div("dash=={}".format(dash.__version__)),
            html.Div("dash-dangerously-set-inner-html=={}".format(dash_dangerously_set_inner_html.__version__)),
            html.Div("mudata=={}".format(mudata.__version__)),
            html.Div("pandas=={}".format(pd.__version__)),
            html.Div(id='package-versions', className='section', style={'width': '75%', 'display': 'inline-block', 'vertical-align': 'top'}),
        ]),
    ])
], fluid=True, className="p-4")


# Callback to populate the covariate dropdown based on the selected dim reduction
@callback(
    Output('covariate-selector', 'options'),
    Input('dim-reduction-selector', 'value')
)
def update_covariate_options(selected_dim_reduction):
    if selected_dim_reduction is None:
        return []
    
    # Get all the covariates from the dataframe that are not prefixed with selected_dim_reduction
    covariate_options = [{'label': col, 'value': col} for col in dim_reduce[selected_dim_reduction].columns if not col.startswith(selected_dim_reduction)]
    return covariate_options

import numpy as np
import matplotlib.pyplot as plt
tab20 = plt.get_cmap('tab20').colors

# Generate a larger categorical colormap using matplotlib
def generate_large_colormap(num_colors):
    cmap = plt.get_cmap('tab20b', num_colors)
    colors = [f'rgba({int(r*255)}, {int(g*255)}, {int(b*255)}, 0.8)' for r, g, b, _ in cmap(np.linspace(0, 1, num_colors))]
    return colors

# Helper function to map categories to colors
def map_categories_to_colors(categories):
    unique_categories = sorted(categories.unique())
    colors = generate_large_colormap(len(unique_categories))
    color_map = {category: color for category, color in zip(unique_categories, colors)}
    return categories.map(color_map).tolist(), color_map

# Callback to update the scatter plot and legend based on selected dim reduction and covariate
@dash.callback(
    [Output('dim-reduction-scatter', 'figure'),
     Output('legend', 'children')],
    [Input('dim-reduction-selector', 'value'), Input('covariate-selector', 'value')]
)
def update_dim_reduction_plot(selected_dim_reduction, selected_covariate):
    if selected_dim_reduction is None or selected_covariate is None:
        return go.Figure(), html.Div()  # Return an empty figure and legend placeholder

    # Extract the relevant data for plotting (they will be selected_dim_reduction_0 and selected_dim_reduction_1)
    data = dim_reduce[selected_dim_reduction][[f'{selected_dim_reduction}_0', f'{selected_dim_reduction}_1', selected_covariate]]
    data.columns = ['X', 'Y', selected_covariate]

    # Map covariate categories to more colors
    colors, color_map = map_categories_to_colors(data[selected_covariate])

    # Plot using the scatterplot function
    fig = scatterplot(
        data=data,
        x_column="X",
        y_column="Y",
        sorted=False,
        title='',
        x_axis_title=f'{selected_dim_reduction}1',
        y_axis_title=f'{selected_dim_reduction}2',
        colors=colors
    )

    # Generate the legend
    legend_items = [html.Div([
                        html.Span(style={'backgroundColor': color, 'display': 'inline-block', 'width': '20px', 'height': '20px', 'marginRight': '10px'}),
                        html.Span(category)
                    ], style={'marginBottom': '5px'}) for category, color in color_map.items()]

    legend = html.Div(legend_items, style={'display': 'flex', 'flexDirection': 'column'})
    
    return fig, legend
