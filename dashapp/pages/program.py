import os
import dash
import pickle
import pandas as pd
import dash_bootstrap_components as dbc
from dash import html, dcc, callback, Input, Output
from dash import dash_table
from plot import scatterplot, barplot, lollipop_plot, boxplot
from dash import DiskcacheManager
import diskcache
import plotly.express as px
import numpy as np
import plotly.graph_objects as go
from dash import Output, Input, State
ANNOTATIONS_FILE = "/cellar/users/aklie/opt/gene_program_evaluation/dashapp/example_data/iPSC_EC_evaluations/annotations.csv"


# Register the page
dash.register_page(__name__, order=3)

# Get the app and cache
app = dash.get_app()
cache = diskcache.Cache("./.cache")
results = cache.get("results")

# Get the first data type and its corresponding second-level keys as default
default_data_type = "explained_variance_ratios"
default_run = list(results[default_data_type].keys())[0]

# Get the first program from the data type and run as default
programs = sorted(list(results[default_data_type][default_run]["program_name"].astype(str).unique()))
default_program = programs[0]

# Get the columns that don't have the selected dim reduction prefix
obsms = results["obsms"]
default_dim_reduction = 'X_umap' if 'X_umap' in obsms else list(obsms.keys())[0]
no_reduce_columns = [col for col in obsms[default_dim_reduction].columns if not col.startswith(default_dim_reduction)]
default_covariate = 'sample' if 'sample' in obsms[default_dim_reduction] else no_reduce_columns[0]
print(f"Default covariate: {default_covariate}")

# Create the layout using tabs for each section
layout = dbc.Container([
    
    # Title
    html.H1("Program Analysis", className="mb-4"),

    # Dropdown to select the run
    dbc.Row([
        dbc.Col([
            html.Label("Select Run"),
            dcc.Dropdown(
                id='run-selector',
                options=[{'label': run, 'value': run} for run in results[default_data_type].keys()],
                value=default_run
            )
        ], width=6),
    ], className="mb-4"),

    # Dropdown to select the program
    dbc.Row([
        dbc.Col([
            html.Label("Select Program"),
            dcc.Dropdown(
                id='program-selector',
                options=[{'label': program, 'value': program} for program in programs],
                value=default_program
            )
        ], width=6),
    ], className="mb-4"),

    # Tabs for different sections of analysis
    dcc.Tabs([

        # Tab 1: Gene Loadings
        dcc.Tab(label='Gene Loadings', children=[
            html.H2("Gene Loadings", className="mt-4 mb-3"),
            html.P("Scatter plot showing top genes by loading for this program. Consider adding links to external gene resources.", className="mb-3"),
            dcc.Input(
                id='gene-loadings-n',
                type='number',
                value=25,
                min=1,
                max=100,
                step=1,
                placeholder="Enter number of genes to display"
            ),
            dcc.Graph(id='gene-loadings-plot'),
        ]),

        # Tab 2: Gene Set Enrichment
        dcc.Tab(label='Gene Set Enrichment', children=[
            html.H2("Gene Set Enrichment Results", className="mt-4 mb-3"),
            
            dbc.Row([
                dbc.Col([
                    html.H3("Table of Enriched Terms"),
                    dcc.Input(
                        id='enriched-terms-threshold',
                        type='number',
                        value=0.05,
                        min=0,
                        max=1,
                        step=0.01,
                        placeholder="Enter significance threshold"
                    ),
                    html.Div(id='enriched-terms-table', className="mb-4"),
                ], width=12),
            ]),

            dbc.Row([
                dbc.Col([
                    html.H3("Volcano Plot of Terms"),
                    dcc.Graph(id='volcano-plot-terms'),
                ], width=12),
            ]),
        ]),

        # Tab 3: Motif Enrichment
        dcc.Tab(label='Motif Enrichment', children=[
            html.H2("Motif Enrichment Results", className="mt-4 mb-3"),
            dbc.Row([
                dbc.Col([
                    html.H3("Table of Enriched Motifs"),
                    dcc.Input(
                        id='enriched-motifs-threshold',
                        type='number',
                        value=0.05,
                        min=0,
                        max=1,
                        step=0.01,
                        placeholder="Enter significance threshold"
                    ),
                    html.Div(id='enriched-motifs-table', className="mb-4"),
                ], width=12),
            ]),

            dbc.Row([
                dbc.Col([
                    html.H3("Volcano Plot of Motifs"),
                    dcc.Graph(id='volcano-plot-motifs'),
                ], width=12),
            ]),
            html.H3("Motif Logo for selected motif", className="mt-3 mb-3"),
            dcc.Graph(id='motif-logo'),
        ]),

        # Tab 4: Trait Enrichment
        dcc.Tab(label='Trait Enrichment', children=[
            html.H2("Trait Enrichment Results", className="mt-4 mb-3"),
            dbc.Row([
                dbc.Col([
                    html.H3("Table of Enriched Trait Terms"),
                    dcc.Input(
                        id='enriched-traits-threshold',
                        type='number',
                        value=0.05,
                        min=0,
                        max=1,
                        step=0.01,
                        placeholder="Enter significance threshold"
                    ),
                    html.Div(id='enriched-traits-table', className="mb-4"),
                ], width=12),
            ]),
            
            dbc.Row([
                dbc.Col([
                    html.H3("Volcano Plot of Traits"),
                    dcc.Graph(id='volcano-plot-traits'),
                ], width=12),
            ]),

            html.H3("Phewas Plot for Selected Program", className="mt-3 mb-3"),
            dcc.Graph(id='phewas-binary-program-plot'),

            html.H3("Phewas Plot for Continuous Traits", className="mt-3 mb-3"),
            dcc.Graph(id='phewas-continuous-program-plot'),
        ]),

        # Tab 5: Categorical Association
        dcc.Tab(label='Categorical Association', children=[
            html.H2("Categorical Association Analysis", className="mt-4 mb-3"),
            
            # Dropdown for selecting covariate
            dbc.Row([
                dbc.Col([
                    html.Label("Select Covariate"),
                    dcc.Dropdown(
                        id='program-covariate-selector',
                        options=["sample"],  # Will be populated based on dim reduction
                        value=default_covariate,
                        clearable=False
                    )
                ], width=6),
            ], className="mb-4"),

            # Boxplot
            dbc.Row([
                dbc.Col([
                    html.H3("Boxplot of Categorical Association"),
                    dcc.Graph(id='categorical-association-program-plot'),
                ], width=12),
            ]),

            # Dropdowns for selecting dim reduction
            dbc.Col([
                    html.Label("Select Dimensionality Reduction"),
                    dcc.Dropdown(
                        id='program-dim-reduction-selector',
                        options=[{'label': key, 'value': key} for key in obsms.keys()],
                        value=default_dim_reduction,
                        clearable=False
                    )
                ], width=6),

            # Scatter plot and legend side by side
            dbc.Row([
                # Plot occupies half the space
                dbc.Col([
                    html.H3("Scatter Plot of Dimensionality Reduction"),
                    dcc.Graph(id='program-dim-reduction-scatter', className="mt-3"),
                ], width=6),  # Half the width
                # Same plot but with selected covariate
                dbc.Col([
                    html.H3("Scatter Plot of Dimensionality Reduction with Covariate"),
                    dcc.Graph(id='program-dim-reduction-covariate-scatter', className="mt-3"),
                ], width=6),  # Half the width
            ]),
        ]),

        # Tab 6: Perturbation
        dcc.Tab(label='Perturbation Association', children=[
            html.H2("Perturbation Association Results", className="mt-4 mb-3"),

            dbc.Row([
                dbc.Col([
                    html.H3("Table of Associated Perturbations"),
                    dcc.Input(
                        id='perturbation-threshold',
                        type='number',
                        value=0.05,
                        min=0,
                        max=1,
                        step=0.01,
                        placeholder="Enter significance threshold"
                    ),
                    html.Div(id='perturbation-table', className="mb-4"),
                ], width=12),
            ]),

            dbc.Row([
                dbc.Col([
                    html.H3("Volcano Plot of Perturbations"),
                    dcc.Graph(id='volcano-plot-perturbations'),
                ], width=12),
            ]),
        ]),

        # Tab 7: Annotations
        dcc.Tab(label='Add annotations', children=[
            html.H2("Type in Annotation", className="mt-4 mb-3"),
            dbc.Row([
                dbc.Col([
                    dcc.Textarea(
                        id='annotation-text',
                        placeholder="Enter annotations or notes here...",
                        style={'width': '100%', 'height': 100},
                    ),
                ], width=12),
            ], className="mb-4"),
            dbc.Row([
                dbc.Col([
                    dbc.Button('Submit Annotation', id='submit-annotation', color='primary', className='mt-3'),
                    html.Div(id='annotation-save-status', className='mt-3'),
                ], width=12),
            ]),
        ]),

    ])
], fluid=True, className="p-4")


# Callback for gene loadings plot
@callback(
    Output('gene-loadings-plot', 'figure'),
    [Input('run-selector', 'value'), Input('program-selector', 'value'), Input('gene-loadings-n', 'value')]
)
def update_gene_loadings_plot(
    selected_run, 
    selected_program, 
    n,
    debug=False,
):
    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program}, n: {n}")

    # Assuming we want to plot something from the selected run
    data_to_plot = results["loadings"][selected_run].loc[selected_program].to_frame(name="loadings").reset_index()

    # Select top n genes by loading
    data_to_plot = data_to_plot.sort_values(by="loadings", ascending=False).head(n).reset_index(drop=True)

    # Define colors based on the sign of the loadings
    marker_colors = ['blue' if val > 0 else 'red' for val in data_to_plot['loadings']]
    line_colors = ['blue' if val > 0 else 'red' for val in data_to_plot['loadings']]

    # Plot the data
    fig = lollipop_plot(
        data=data_to_plot,
        x_column='gene_name',
        y_column='loadings',
        title='Gene Loadings for Selected Program',
        x_axis_title='Gene',
        y_axis_title='Loadings',
        marker_colors=marker_colors,
        line_colors=line_colors,
    )

    return fig


# Callback for geneset enrichment table
@callback(
    Output('enriched-terms-table', 'children'),
    [Input('run-selector', 'value'), Input('program-selector', 'value'), Input('enriched-terms-threshold', 'value')]
)
def update_enriched_terms_table(
    selected_run, 
    selected_program, 
    sig_threshold,
    debug=False,
    categorical_var = "program_name",
    sig_var = "adj_pval",
):  

    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program}, sig_threshold: {sig_threshold}")

    # Retrieve the relevant data
    data_to_plot = results['geneset_enrichments'].get(selected_run, pd.DataFrame())

    # Make sure x-axis is string
    data_to_plot[categorical_var] = data_to_plot[categorical_var].astype(str)

    # Filter data for the selected program
    data_to_plot = data_to_plot[data_to_plot[categorical_var] == selected_program]

    # Filter data for significant terms
    data_to_plot = data_to_plot[data_to_plot[sig_var] < sig_threshold]

    # Sort data by adj_pval
    data_to_plot = data_to_plot.sort_values(by=sig_var)

    if data_to_plot.empty:
        print("No significant terms found.")
        return html.Div("No significant enriched terms found.")

    # Create a DataTable
    table = dash_table.DataTable(
        id='enriched-terms-data-table',
        columns=[{"name": col, "id": col} for col in data_to_plot.columns],
        data=data_to_plot.to_dict('records'),
        filter_action="native",  # Enables filtering
        sort_action="native",    # Enables sorting
        page_size=10,            # Number of rows per page
        style_table={'overflowX': 'auto'},  # Scrollable horizontally if too wide
        style_cell={
            'textAlign': 'left',  # Align text to the left
            'padding': '5px',
        },
        style_header={
            'backgroundColor': 'rgb(230, 230, 230)',
            'fontWeight': 'bold'
        },
    )

    return table


# Callback for volcano plot of terms
@callback(
    Output('volcano-plot-terms', 'figure'),
    [Input('run-selector', 'value'), Input('program-selector', 'value')]
)
def update_volcano_plot_terms(
    selected_run,
    selected_program,
    categorical_var = "term",
    sig_var = "adj_pval",
    debug=False,
):

    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program}")

    # Assuming we want to plot something from the selected run
    data_to_plot = results['geneset_enrichments'][selected_run]

    # Filter data for the selected program
    data_to_plot = data_to_plot.query(f"program_name == '{selected_program}'")

    # -log10 transform adj_pvals
    data_to_plot[f'-log10({sig_var})'] = -np.log10(data_to_plot[sig_var])
    
    # Make sure x-axis is string
    data_to_plot[categorical_var] = data_to_plot[categorical_var].astype(str)

    # Example plot using Plotly (replace with your actual plotting code)
    fig = scatterplot(data=data_to_plot,
        x_column=categorical_var,
        y_column=f'-log10({sig_var})',
        title='',
        x_axis_title=categorical_var,
        y_axis_title=f'-log10({sig_var})',
    )
    return fig


# Callback for motif enrichment table
@callback(
    Output('enriched-motifs-table', 'children'),
    [Input('run-selector', 'value'), Input('program-selector', 'value'), Input('enriched-motifs-threshold', 'value')]
)
def update_enriched_motifs_table(
    selected_run, 
    selected_program, 
    sig_threshold,
    debug=False,
    categorical_var = "program_name",
    sig_var = "pval",
):
  
    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program}, sig_threshold: {sig_threshold}")

    # Retrieve the relevant data
    data_to_plot = results['motif_enrichments'].get(selected_run, pd.DataFrame())

    # Make sure x-axis is string
    data_to_plot[categorical_var] = data_to_plot[categorical_var].astype(str)

    # Filter data for the selected program
    data_to_plot = data_to_plot[data_to_plot[categorical_var] == selected_program]

    # Filter data for significant terms
    data_to_plot = data_to_plot[data_to_plot[sig_var] < sig_threshold]
    
    # Sort data by adj_pval
    data_to_plot = data_to_plot.sort_values(by=sig_var)

    if data_to_plot.empty:
        print("No significant motifs found.")
        return html.Div("No significant enriched motifs found.")

    # Create a DataTable
    table = dash_table.DataTable(
        id='enriched-motifs-data-table',
        columns=[{"name": col, "id": col} for col in data_to_plot.columns],
        data=data_to_plot.to_dict('records'),
        filter_action="native",  # Enables filtering
        sort_action="native",    # Enables sorting
        page_size=10,            # Number of rows per page
        style_table={'overflowX': 'auto'},  # Scrollable horizontally if too wide
        style_cell={
            'textAlign': 'left',  # Align text to the left
            'padding': '5px',
        },
        style_header={
            'backgroundColor': 'rgb(230, 230, 230)',
            'fontWeight': 'bold'
        },
    )

    return table


# Callback for volcano plot of motifs
@callback(
    Output('volcano-plot-motifs', 'figure'),
    [Input('run-selector', 'value'), Input('program-selector', 'value')]
)
def update_volcano_plot_motifs(
    selected_run,
    selected_program,
    categorical_var = "motif",
    sig_var = "pval",
    debug=False,
):

    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program}")

    # Assuming we want to plot something from the selected run
    data_to_plot = results['motif_enrichments'][selected_run]

    # Filter data for the selected program
    data_to_plot = data_to_plot.query(f"program_name == '{selected_program}'")

    # -log10 transform adj_pvals
    data_to_plot[f'-log10({sig_var})'] = -np.log10(data_to_plot[sig_var])

    # Make sure x-axis is string
    data_to_plot[categorical_var] = data_to_plot[categorical_var].astype(str)

    # Example plot using Plotly (replace with your actual plotting code)
    fig = scatterplot(
        data=data_to_plot,
        x_column=categorical_var,
        y_column=f'-log10({sig_var})',
        title='',
        x_axis_title=categorical_var,
        y_axis_title=f'-log10({sig_var})',
    )
    return fig


# Callback for trait enrichment table
@callback(
    Output('enriched-traits-table', 'children'),
    [Input('run-selector', 'value'), Input('program-selector', 'value'), Input('enriched-traits-threshold', 'value')]
)
def update_enriched_traits_table(
    selected_run, 
    selected_program,
    sig_threshold,
    categorical_var = "program_name",
    sig_var = "adj_pval",
    debug=False,
):

    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program}, sig_threshold: {sig_threshold}")

    # Retrieve the relevant data
    data_to_plot = results['trait_enrichments'].get(selected_run, pd.DataFrame())

    # Make sure x-axis is string
    data_to_plot[categorical_var] = data_to_plot[categorical_var].astype(str)

    # Filter data for the selected program
    data_to_plot = data_to_plot[data_to_plot[categorical_var] == selected_program]

    # Filter data for significant terms
    data_to_plot = data_to_plot[data_to_plot[sig_var] < sig_threshold]

    # Sort data by adj_pval
    data_to_plot = data_to_plot.sort_values(by=sig_var)

    if data_to_plot.empty:
        print("No significant traits found.")
        return html.Div("No significant enriched traits found.")

    # Create a DataTable
    table = dash_table.DataTable(
        id='enriched-traits-data-table',
        columns=[{"name": col, "id": col} for col in data_to_plot.columns],
        data=data_to_plot.to_dict('records'),
        filter_action="native",  # Enables filtering
        sort_action="native",    # Enables sorting
        page_size=10,            # Number of rows per page
        style_table={'overflowX': 'auto'},  # Scrollable horizontally if too wide
        style_cell={
            'textAlign': 'left',  # Align text to the left
            'padding': '5px',
        },
        style_header={
            'backgroundColor': 'rgb(230, 230, 230)',
            'fontWeight': 'bold'
        },
    )

    return table


# Callback for volcano plot of traits
@callback(
    Output('volcano-plot-traits', 'figure'),
    [Input('run-selector', 'value'), Input('program-selector', 'value')]
)
def update_volcano_plot_traits(
    selected_run,
    selected_program,
    categorical_var = "trait_reported",
    sig_var = "adj_pval",
    debug=False,
):

    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program}")

    # Assuming we want to plot something from the selected run
    data_to_plot = results['trait_enrichments'][selected_run]

    # Filter data for the selected program
    data_to_plot = data_to_plot.query(f"program_name == '{selected_program}'")

    # -log10 transform adj_pvals
    data_to_plot[f'-log10({sig_var})'] = -np.log10(data_to_plot[sig_var])

    # Make sure x-axis is string
    data_to_plot[categorical_var] = data_to_plot[categorical_var].astype(str)

    # Example plot using Plotly (replace with your actual plotting code)
    fig = scatterplot(
        data=data_to_plot,
        x_column=categorical_var,
        y_column=f'-log10({sig_var})',
        title='',
        x_axis_title=categorical_var,
        y_axis_title=f'-log10({sig_var})',
    )
    return fig


# Callback for phewas binary plot
@callback(
    Output('phewas-binary-program-plot', 'figure'),
    [Input('run-selector', 'value'), Input('program-selector', 'value')]
)
def update_phewas_binary_plot(
    selected_run,
    selected_program,
    debug=False,
):

    # Assuming we want to plot something from the selected run
    data_to_plot = results['trait_enrichments'][selected_run]

    # Filter data for the selected program
    data_to_plot = data_to_plot.query(f"program_name == '{selected_program}'")

    # Filter data for binary traits
    data_to_plot = data_to_plot.query("trait_category != 'measurement'")

    fig = px.scatter(
        data_to_plot.query("trait_category != 'measurement'"),
        x='trait_reported',
        y='-log10(adj_pval)',
        color='trait_category',
        title="",
        hover_data=["program_name", "trait_reported", "trait_category", "adj_pval", "genes", "study_id", "pmid"]
    )

    # Customize layout
    fig.update_layout(
        xaxis_title='trait_reported',
        yaxis_title='-log10(adj_pval)',
        yaxis=dict(tickformat=".1f"),
        width=1000,
        height=800,
        xaxis_tickfont=dict(size=4),
        plot_bgcolor='rgba(0,0,0,0)',  # Transparent background
    )

    # Add horizontal dashed line for significance threshold
    fig.add_hline(
        y=-np.log10(0.05), 
        line_dash="dash",
        annotation_text='Significance Threshold (0.05)', 
        annotation_position="top right"
    )

    return fig


# Callback for phewas continuous plot
@callback(
    Output('phewas-continuous-program-plot', 'figure'),
    [Input('run-selector', 'value'), Input('program-selector', 'value')]
)
def update_phewas_continuous_plot(
    selected_run,
    selected_program,
    debug=False,
):
    
        # Assuming we want to plot something from the selected run
        data_to_plot = results['trait_enrichments'][selected_run]
    
        # Filter data for the selected program
        data_to_plot = data_to_plot.query(f"program_name == '{selected_program}'")
    
        # Filter data for continuous traits
        data_to_plot = data_to_plot.query("trait_category == 'measurement'")
    
        fig = px.scatter(
            data_to_plot,
            x='trait_reported',
            y='-log10(adj_pval)',
            color='trait_category',
            title="",
            hover_data=["program_name", "trait_reported", "trait_category", "adj_pval", "genes", "study_id", "pmid"]
        )
    
        # Customize layout
        fig.update_layout(
            xaxis_title='trait_reported',
            yaxis_title='-log10(adj_pval)',
            yaxis=dict(tickformat=".1f"),
            width=1000,
            height=800,
            xaxis_tickfont=dict(size=4),
            plot_bgcolor='rgba(0,0,0,0)',  # Transparent background
        )
    
        # Add horizontal dashed line for significance threshold
        fig.add_hline(
            y=-np.log10(0.05), 
            line_dash="dash",
            annotation_text='Significance Threshold (0.05)', 
            annotation_position="top right"
        )
    
        return fig


# Callback for covariate association plot
@callback(
    Output('categorical-association-program-plot', 'figure'),
    [Input('run-selector', 'value'), Input('program-selector', 'value'), Input('program-covariate-selector', 'value')]
)
def update_covariate_association_plot(
    selected_run, 
    selected_program, 
    selected_covariate, 
    debug=False
):

    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program}, covariate: {selected_covariate}")

    if selected_covariate is None:
        return go.Figure()
    
    curr_obs = results['obs'][selected_run]
    curr_obs_membership = results['obs_memberships'][selected_run]
    concat_df = pd.concat([curr_obs_membership[selected_program], curr_obs[selected_covariate]], axis=1)
    
    fig = boxplot(
        concat_df, 
        x_column=selected_program, 
        y_column=selected_covariate, 
        title="",
        x_axis_title="Cell membership value",
        y_axis_title=selected_covariate,
    )

    return fig


# Callback to update the scatter plot and legend based on selected dim reduction and covariate
@callback(
    Output('program-dim-reduction-scatter', 'figure'),
    [Input('program-dim-reduction-selector', 'value'), Input('run-selector', 'value'), Input('program-selector', 'value')]
)
def update_program_dim_reduction_plot(
    selected_dim_reduction,
    selected_run,
    selected_program,
    debug=False
):
    
    if debug:
        print(f"Selected dim reduction: {selected_dim_reduction}, selected run: {selected_run}, selected program: {selected_program}")

    if selected_dim_reduction is None:
        return go.Figure(), html.Div()  # Return an empty figure and legend placeholder

    # Get vector of continuous values for the program membership
    curr_obs = results['obs'][selected_run]
    print(curr_obs)
    curr_obs_membership = results['obs_memberships'][selected_run][selected_program]
    print(f"Current obs membership: {curr_obs_membership}")

    # Extract the relevant data for plotting (they will be selected_dim_reduction_0 and selected_dim_reduction_1)
    data = obsms[selected_dim_reduction][[f'{selected_dim_reduction}_0', f'{selected_dim_reduction}_1']]
    data.columns = ['X', 'Y']

    # Plot using the scatterplot function
    fig = scatterplot(
        data=data,
        x_column="X",
        y_column="Y",
        sorted=False,
        title='',
        x_axis_title=f'{selected_dim_reduction}1',
        y_axis_title=f'{selected_dim_reduction}2',
        colors=curr_obs_membership
    )
    
    return fig


# Callback to update the scatter plot and legend based on selected dim reduction and covariate
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

@callback(
    Output('program-dim-reduction-covariate-scatter', 'figure'),
    [Input('program-dim-reduction-selector', 'value'), Input('run-selector', 'value'), Input('program-covariate-selector', 'value')]
)
def update_program_dim_reduction_covariate_plot(
    selected_dim_reduction,
    selected_run,
    selected_covariate,
    debug=True
):
    
    if debug:
        print(f"Selected dim reduction: {selected_dim_reduction}, selected run: {selected_run}, selected covariate: {selected_covariate}")

    if selected_dim_reduction is None or selected_covariate is None:
        return go.Figure(), html.Div()
    
    # Extract the relevant data for plotting (they will be selected_dim_reduction_0 and selected_dim_reduction_1)
    data = obsms[selected_dim_reduction][[f'{selected_dim_reduction}_0', f'{selected_dim_reduction}_1', selected_covariate]]
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

    return fig


# Callback for perturbation table
@callback(
    Output('perturbation-table', 'children'),
    [Input('run-selector', 'value'), Input('program-selector', 'value'), Input('perturbation-threshold', 'value')]
)
def update_perturbation_table(
    selected_run, 
    selected_program, 
    sig_threshold,
    debug=True,
    categorical_var = "target_name",
    sig_var = "pval",
):
    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program}, sig_threshold: {sig_threshold}")

    # Retrieve the relevant data
    data_to_plot = results['perturbation_associations'].get(selected_run, pd.DataFrame())

    # Make sure x-axis is string
    data_to_plot[categorical_var] = data_to_plot[categorical_var].astype(str)

    # Filter data for the selected program
    data_to_plot = data_to_plot[data_to_plot['program_name'].astype(str) == selected_program]

    # Filter data for significant terms
    data_to_plot = data_to_plot[data_to_plot[sig_var] < sig_threshold]
    
    # Sort data by adj_pval
    data_to_plot = data_to_plot.sort_values(by=sig_var)

    if data_to_plot.empty:
        print("No significant perturbation associations found.")
        return html.Div("No significant perturbation associations found.")

    # Create a DataTable
    table = dash_table.DataTable(
        id='perturbation-data-table',
        columns=[{"name": col, "id": col} for col in data_to_plot.columns],
        data=data_to_plot.to_dict('records'),
        filter_action="native",  # Enables filtering
        sort_action="native",    # Enables sorting
        page_size=10,            # Number of rows per page
        style_table={'overflowX': 'auto'},  # Scrollable horizontally if too wide
        style_cell={
            'textAlign': 'left',  # Align text to the left
            'padding': '5px',
        },
        style_header={
            'backgroundColor': 'rgb(230, 230, 230)',
            'fontWeight': 'bold'
        },
    )

    return table


# Callback for perturbation association plot
@callback(
    Output('volcano-plot-perturbations', 'figure'),
    [Input('run-selector', 'value'), Input('program-selector', 'value')]
)
def update_perturbation_association_plot(
    selected_run, 
    selected_program,
    categorical_var = "target_name",
    sig_var = "pval",
    sig_threshold = 0.05,
    debug=True
):
    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program}")

    # Retrieve the relevant data
    data_to_plot = results['perturbation_associations'].get(selected_run, pd.DataFrame())

    # Make sure x-axis is string
    data_to_plot[categorical_var] = data_to_plot[categorical_var].astype(str)

    # Filter data for the selected program
    data_to_plot = data_to_plot[data_to_plot['program_name'].astype(str) == selected_program]

    # -log10 transform adj_pvals
    data_to_plot[f'-log10({sig_var})'] = -np.log10(data_to_plot[sig_var])

    # Filter data for significant terms
    data_to_plot = data_to_plot[data_to_plot[sig_var] < sig_threshold]

    # Sort data by adj_pval
    data_to_plot = data_to_plot.sort_values(by=sig_var)

    # Example plot using Plotly (replace with your actual plotting code)
    fig = scatterplot(
        data=data_to_plot,
        x_column=categorical_var,
        y_column=f'-log10({sig_var})',
        title='',
        x_axis_title=categorical_var,
        y_axis_title=f'-log10({sig_var})',
    )

    return fig


@app.callback(
    Output('annotation-save-status', 'children'),
    Input('submit-annotation', 'n_clicks'),
    State('program-selector', 'value'),
    State('annotation-text', 'value')
)
def save_annotation(n_clicks, selected_program, annotation_text):
    if n_clicks is None:
        return ""
    
    # Read the existing annotations file
    if os.path.exists(ANNOTATIONS_FILE):
        df = pd.read_csv(ANNOTATIONS_FILE)
    else:
        df = pd.DataFrame(columns=['Program', 'Annotation'])
    
    # Update or append the annotation
    if selected_program in df['Program'].values:
        df.loc[df['Program'] == selected_program, 'Annotation'] = annotation_text
    else:
        df = pd.concat([df, pd.DataFrame({'Program': [selected_program], 'Annotation': [annotation_text]})])
    
    # Save the updated DataFrame back to the CSV file
    df.to_csv(ANNOTATIONS_FILE, index=False)
    
    return f"Annotation for '{selected_program}' has been saved."