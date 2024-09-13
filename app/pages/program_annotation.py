import os
import dash
import pickle
import pandas as pd
import dash_bootstrap_components as dbc
from dash import html, dcc, callback, Input, Output
from dash import dash_table
from plot import scatterplot, volcano_plot, lollipop_plot, boxplot, fig_to_uri, scatterplot_static
from utils import map_categories_to_colors
import diskcache
import plotly.express as px
import numpy as np
import plotly.graph_objects as go
from dash import Output, Input, State
import datetime
import matplotlib.pyplot as plt


# Register the page
dash.register_page(__name__, order=3)

# Get the app and cache
app = dash.get_app()
cache = diskcache.Cache("./.cache")
results = cache.get("results")

# Get
path_report_out = results["path_report_out"]
path_mdata = results["path_mdata"]
path_evaluation_outs = results["path_evaluation_outs"]
data_key = results["data_key"]
annotations_loc = results["annotations_loc"]

# Get the first data type and its corresponding second-level keys as default
default_data_type = "explained_variance_ratios"
default_run = list(results[default_data_type].keys())[0]

# Get the first program from the data type and run as default
# if the programs are numbers sort them in numerical order
programs = sorted(list(results[default_data_type][default_run]["program_name"].astype(str).unique()))
if programs[0].isdigit():
    programs = sorted(programs, key=int)
default_program = programs[0]

# Grab the dimensionality reduction data and categorical keys
obsms = results["obsms"]
default_obsm = 'X_umap' if 'X_umap' in obsms else list(obsms.keys())[0]

# Get the columns that don't have the selected dim reduction prefix
categorical_keys = results["categorical_keys"]
continuous_keys = results["continuous_keys"]
covariate_keys = categorical_keys + continuous_keys
default_covariate = categorical_keys[0] if categorical_keys else None

# Grab the groupings
perturbation_association_stratification_key = results['perturbation_association_stratification_key']
motif_enrichment_stratification_key = results['motif_enrichment_stratification_key']

# Get all the possible trait_category values for trait enrichment
trait_categories = results['trait_enrichments'][default_run]['results'][list(results['trait_enrichments'][default_run]['results'].keys())[0]]['trait_category'].unique()
trait_categories_to_show = [trait_category for trait_category in trait_categories if trait_category not in ["biological process", "disease of ear", "phenotype", "injury, poisoning or other complication"]]

# Create the layout using tabs for each section
layout = dbc.Container([
    
    # Title
    html.H1("Program Analysis", className="mb-4"),
    html.P("This page is designed for analysis and annotation of individual programs. "
           "It provide several layers of information about the program, "
           "including the genes at the top of the loadings, the enriched gene sets, "
           "the enriched motifs, the enriched traits, the perturbation associations, "
           "and the categorical associations, which can aid in the annotation of the program. "
           "The 'Add Annotations' tab allows you to add notes or comments to the program. "
           "That are saved and can be viewed later.", className="mb-4"),

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
            html.H2("Gene loadings for selected program", className="mt-4 mb-3"),
            html.P("This tab provides a look at the top gene loadings for the selected program. "
                   "Use the input box to specify the number of genes you want to display.", className="mb-4"),

            # Input for number of genes to display
            dcc.Input(
                id='gene-loadings-n',
                type='number',
                value=25,
                min=1,
                max=100,
                step=1,
                placeholder="Enter number of genes to display"
            ),

            # gene loadings plot
            dcc.Graph(id='gene-loadings-plot'),
            
            # gene loadings table
            html.Div(id='gene-loadings-table', className="mb-4"),
        ]),

        # Tab 2: Gene Set Enrichment
        dcc.Tab(label='Gene Set Enrichment', children=[
            html.H2("Enriched terms for selected program", className="mt-4 mb-3"),
            html.P("This tab provides a look at the enriched gene sets for the selected program. "
                   "Use the input box to specify the significance threshold you want to use.", className="mb-4"),
            
            dbc.Row([
                
                # Drop down for library
                dbc.Col([
                    html.Label("Select Library"),
                    dcc.Dropdown(
                        id='table-genesets-library-selector',
                        options=[{'label': library, 'value': library} for library in list(set(results['geneset_enrichments'][default_run]['libraries']))],
                        value=list(set(results['geneset_enrichments'][default_run]['libraries']))[0]
                    )
                ], width=6),

                # Drop down for method
                dbc.Col([
                    html.Label("Select Method"),
                    dcc.Dropdown(
                        id='table-genesets-method-selector',
                        options=[{'label': method, 'value': method} for method in list(set(results['geneset_enrichments'][default_run]['methods']))],
                        value=list(set(results['geneset_enrichments'][default_run]['methods']))[0]
                    )
                ], width=6),

            ]),

            dbc.Row([

                # Pvalue input
                dbc.Col([
                    html.Label("Select Significance Threshold"),
                    dcc.Input(
                        id='enriched-terms-threshold',
                        type='number',
                        value=0.05,
                        min=0,
                        max=1,
                        placeholder="Enter significance threshold"
                    ),
                ], width=6),

                # Enrichment input
                dbc.Col([
                    html.Label("Select Enrichment Threshold"),
                    dcc.Input(
                        id='enriched-terms-enrichment-threshold',
                        type='number',
                        value=0,
                        min=0,
                        max=1e6,
                        placeholder="Enter enrichment threshold"
                    ),
                ], width=6),

                
            ]),

            # Volcano plot of terms
            dbc.Row([
                dbc.Col([
                    html.H3("Volcano plot of gene set enrichment"),
                    html.P("This plot depicts statistical significance of gene sets in the program. "
                           "on the y-axis and the effect size on the x-axis. The effect size calculation "
                           "is dependent on the method used for the enrichment analysis. If GSEA was used, "
                           "the effect size is the normalized enrichment score (NES). If a Fisher's exact test "
                            "was used, the effect size is the odds ratio."),
                    dcc.Graph(id='volcano-plot-terms'),
                ], width=12),
            ]),

            # Enriched terms table
            dbc.Row([
                dbc.Col([
                    html.H3("Enriched terms"),
                    html.Div(id='enriched-terms-table', className="mb-4"),
                ], width=12),
            ]),

            
        ]),

        # Tab 3: Motif Enrichment
        dcc.Tab(label='Motif Enrichment', children=[
            html.H2("Enriched motifs in promoters/enhancers of genes in selected program", className="mt-4 mb-3"),
            html.P("This tab provides a look at the enriched motifs for the selected program. "
                   "Use the input box to specify the significance threshold you want to use.", className="mb-4"),

            # Threshold input and E/P type dropdown and database dropdown
            dbc.Row([

                # E/P type dropdown
                dbc.Col([
                    html.Label("Select E/P type"),
                    dcc.Dropdown(
                        id='table-motifs-E_P_type-selector',
                        options=[{'label': E_P_types, 'value': E_P_types} for E_P_types in list(set(results['motif_enrichments'][default_run]['E_P_types']))],
                        value=list(set(results['motif_enrichments'][default_run]['E_P_types']))[0]
                    )
                ], width=3),

                # Database dropdown
                dbc.Col([
                    html.Label("Select Motif Database"),
                    dcc.Dropdown(
                        id='table-motifs-database-selector',
                        options=[{'label': databases, 'value': databases} for databases in list(set(results['motif_enrichments'][default_run]['databases']))],
                        value=list(set(results['motif_enrichments'][default_run]['databases']))[0]
                    )
                ], width=3),

                # Test type dropdown
                dbc.Col([
                    html.Label("Select Motif Enrichment Test Type"),
                    dcc.Dropdown(
                        id='table-motifs-test_type-selector',
                        options=[{'label': test_type, 'value': test_type} for test_type in list(set(results['motif_enrichments'][default_run]['test_types']))],
                        value=list(set(results['motif_enrichments'][default_run]['test_types']))[0]
                    )
                ], width=3),

                # Level key dropdown
                dbc.Col([
                    html.Label("Select Level Key"),
                    dcc.Dropdown(
                        id='table-motifs-level_key-selector',
                        options=[{'label': level_key, 'value': level_key} for level_key in sorted(list(set(results['motif_enrichments'][default_run]['level_keys'])))],
                        value=sorted(list(set(results['motif_enrichments'][default_run]['level_keys'])))[0]
                    )
                ], width=3),

            ]),

            # Test type and level key dropdowns
            dbc.Row([

                # Pvalue input
                dbc.Col([
                    html.Label("Select adjusted p-value significance threshold"),
                    dcc.Input(
                        id='enriched-motifs-threshold',
                        type='number',
                        value=0.05,
                        min=0,
                        max=1,
                        placeholder="Enter significance threshold"
                    ),
                ], width=6),

                # Enrichment value input
                dbc.Col([
                    html.Label("Select Enrichment Threshold"),
                    dcc.Input(
                        id='enriched-motifs-enrichment-threshold',
                        type='number',
                        value=0,
                        min=0,
                        max=1e6,
                        placeholder="Enter enrichment threshold"
                    ),
                ], width=6),
                
            ]),

            # Volcano plot of motifs
            dbc.Row([
                dbc.Col([
                    html.H3("Volcano Plot of Motifs"),
                    dcc.Graph(id='volcano-plot-motifs'),
                ], width=12),
            ]),

            # Enriched motifs table
            dbc.Row([
                dbc.Col([
                    html.H3("Enriched motifs"),
                    html.Div(id='enriched-motifs-table', className="mb-4"),
                ], width=12),
            ]),

            # Motif logo
            dbc.Row([
                dbc.Col([
                    html.H3("Motif Logo for selected motif", className="mt-3 mb-3"),
                    dcc.Graph(id='motif-logo'),
                ], width=12),
            ]),
        ]),

        # Tab 4: Trait Enrichment
        dcc.Tab(label='Trait Enrichment', children=[
            
            dbc.Row([
                 
                 # Drop down for library
                dbc.Col([
                    html.Label("Select Database"),
                    dcc.Dropdown(
                        id='table-traits-library-selector',
                        options=[{'label': library, 'value': library} for library in list(set(results['trait_enrichments'][default_run]['databases']))],
                        value=list(set(results['trait_enrichments'][default_run]['databases']))[0]
                    )
                ], width=4),

                # drop down for method
                dbc.Col([
                    html.Label("Select Method"),
                    dcc.Dropdown(
                        id='table-traits-method-selector',
                        options=[{'label': method, 'value': method} for method in list(set(results['trait_enrichments'][default_run]['methods']))],
                        value=list(set(results['trait_enrichments'][default_run]['methods']))[0]
                    )
                ], width=4),

            ]),
            
            # Threshold inputs
            dbc.Row([

                # Pvalue input
                dbc.Col([
                    html.Label("Select Significance Threshold"),
                    dcc.Input(
                        id='enriched-traits-threshold',
                        type='number',
                        value=0.05,
                        min=0,
                        max=1,
                        placeholder="Enter significance threshold"
                    ),
                ], width=6),

                # Enrichment input
                dbc.Col([
                    html.Label("Select Enrichment Threshold"),
                    dcc.Input(
                        id='enriched-traits-enrichment-threshold',
                        type='number',
                        value=0,
                        min=0,
                        max=1e6,
                        placeholder="Enter enrichment threshold"
                    ),
                ], width=6),

            ]),
            
            # Toggles for trait type and trait categories
            dbc.Row([
                
                # Toggle for binary or continuous
                dbc.Col([
                    html.Label("Select Trait Type"),
                    dcc.RadioItems(
                        id='trait-type-selector',
                        options=[
                            {'label': 'Binary', 'value': 'binary'},
                            {'label': 'Continuous', 'value': 'continuous'},
                        ],
                        value='binary',
                        inline=True,
                    ),
                ], width=2),

                # Checkboxes for all possible trait categories
                dbc.Col([
                    html.Label("Select Trait Categories to Display"),
                    dcc.Checklist(
                        id='trait-category-selector',
                        options=[{'label': trait_category, 'value': trait_category} for trait_category in trait_categories],
                        value=trait_categories_to_show,
                        inline=True,
                    ),
                ], width=10),

            ]),

            # Binary or continuous trait enrichment PheWAS plot
            dbc.Row([
                
                # Plot
                dbc.Col([
                    html.H3("PheWAS plot for selected program", className="mb-4"),
                    dcc.Graph(id='phewas-program-plot'),
                ], width=12),

            ]),
            
            # Enriched traits table
            dbc.Row([
                
                dbc.Col([
                    html.H3("Enriched traits"),
                    html.Div(id='enriched-traits-table', className="mb-1"),
                ], width=12),

            ]),

        ]),

        # Tab 5: Categorical Association
        dcc.Tab(label='Categorical Association', children=[
            html.H2("Categorical Association Analysis", className="mt-4 mb-3"),

            # Covariate and obsm dropdowns at the top
            dbc.Row([

                # Dropdown for selecting covariate
                dbc.Col([
                    html.Label("Select Covariate"),
                    dcc.Dropdown(
                        id='program-covariate-selector',
                        options=[{'label': key, 'value': key} for key in categorical_keys],
                        value=default_covariate,
                        clearable=False
                    )
                ], width=6),

                # Dropdown for selecting dim reduction
                dbc.Col([
                    html.Label("Select Dimensionality Reduction"),
                    dcc.Dropdown(
                        id='program-dim-reduction-selector',
                        options=[{'label': key, 'value': key} for key in obsms.keys()],
                        value=default_obsm,
                        clearable=False
                    )
                ], width=6),

            ], className="mb-4"),        

            # Boxplot and scatter plot
            dbc.Row([
                
                # Boxplot of categorical association
                dbc.Col([
                    html.H3("Boxplot of Categorical Association"),
                    dcc.Graph(id='categorical-association-program-plot'),
                ], width=6),

                # Scatter plot of dim reduction
                dbc.Col([
                    html.Div([html.Img(id='program-dim-reduction-scatter')], className="mt-3"),
                ], width=6),  # Half the width
            ]),
        ]),

        # Tab 6: Perturbation
        dcc.Tab(label='Perturbation Association', children=[
            html.H2("Perturbation Association Results", className="mt-4 mb-3"),

            dbc.Row([
                dbc.Col([
                    dcc.Input(
                        id='perturbation-threshold',
                        type='number',
                        value=0.05,
                        min=0,
                        max=1,
                        placeholder="Enter significance threshold"
                    ),
                    html.Div(id='perturbation-table', className="mb-4"),
                ], width=4),


                dbc.Col([
                    html.Label("Select Gene Guide"),
                    dcc.Dropdown(
                        id='table-perturbations-gene_guide-selector',
                        options=[{'label': gene_guide, 'value': gene_guide} for gene_guide in list(set(results['perturbation_associations'][default_run]['gene_guides']))],
                        value=list(set(results['perturbation_associations'][default_run]['gene_guides']))[0]
                    )
                ], width=4),

                dbc.Col([
                    html.Label("Select Level Key"),
                    dcc.Dropdown(
                        id='table-perturbations-level_key-selector',
                        options=[{'label': level_key, 'value': level_key} for level_key in sorted(list(set(results['perturbation_associations'][default_run]['level_keys'])))],
                        value=sorted(list(set(results['perturbation_associations'][default_run]['level_keys'])))[0]
                    )
                ], width=4),
            ]),

            # Table of perturbations
            dbc.Row([
                dbc.Col([
                    html.H3("Perturbation Associations"),
                    html.Div(id='perturbation-table', className="mb-4"),
                ], width=12),
            ]),

            # Volcano plot of perturbations
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
    [
        Input('run-selector', 'value'), 
        Input('program-selector', 'value'), 
        Input('gene-loadings-n', 'value')
    ]
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
        title='',
        x_axis_title='Gene',
        y_axis_title='Loadings',
        marker_colors=marker_colors,
        line_colors=line_colors,
    )

    return fig


# Callback for gene loadings table
@callback(
    Output('gene-loadings-table', 'children'),
    [
        Input('run-selector', 'value'), 
        Input('program-selector', 'value'), 
        Input('gene-loadings-n', 'value')
    ]
)
def update_gene_loadings_table(
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

    # Create a DataTable
    table = dash_table.DataTable(
        id='gene-loadings-data-table',
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


# Callback for geneset enrichment table
@callback(
    Output('enriched-terms-table', 'children'),
    [
        Input('run-selector', 'value'), 
        Input('program-selector', 'value'), 
        Input('enriched-terms-threshold', 'value'),
        Input('enriched-terms-enrichment-threshold', 'value'),
        Input('table-genesets-library-selector', 'value'),
        Input('table-genesets-method-selector', 'value')
    ]
)
def update_enriched_terms_table(
    selected_run, 
    selected_program, 
    sig_threshold,
    enrichment_threshold,
    library,
    method,
    categorical_var = "program_name",
    sig_var = "adj_pval",
    debug=False,
):  

    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program}, sig_threshold: {sig_threshold}")
        print(f"Library: {library}, Method: {method}")

    # Assuming we want to plot something from the selected run
    data_to_plot = results['geneset_enrichments'][selected_run]["results"][f"{library}_{method}"]

    # Make sure x-axis is string
    data_to_plot[categorical_var] = data_to_plot[categorical_var].astype(str)

    # Filter data for the selected program
    data_to_plot = data_to_plot[data_to_plot[categorical_var] == selected_program]

    # Filter data for significant terms
    data_to_plot = data_to_plot[(data_to_plot[sig_var] < sig_threshold) & (data_to_plot["enrichment"] > enrichment_threshold)]

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
    [
        Input('run-selector', 'value'), 
        Input('program-selector', 'value'),
        Input('table-genesets-library-selector', 'value'),
        Input('table-genesets-method-selector', 'value'),
        Input('enriched-terms-threshold', 'value'),
        Input('enriched-terms-enrichment-threshold', 'value')
    ]
)
def update_volcano_plot_terms(
    selected_run,
    selected_program,
    library,
    method,
    sig_threshold,
    enrichment_threshold,
    categorical_var = "term",
    sig_var = "adj_pval",
    enrichment_var = "enrichment",
    hover_data=["program_name", "term", "adj_pval"],
    debug=False,
):

    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program}")
        print(f"Library: {library}, Method: {method}")

    # Assuming we want to plot something from the selected run
    data_to_plot = results['geneset_enrichments'][selected_run]["results"][f"{library}_{method}"]

    # Filter data for the selected program
    data_to_plot = data_to_plot.query(f"program_name == '{selected_program}'")

    # -log10 transform adj_pvals
    data_to_plot[f'-log10({sig_var})'] = -np.log10(data_to_plot[sig_var])
    
    # Make sure x-axis is string
    data_to_plot[categorical_var] = data_to_plot[categorical_var].astype(str)

    # -log10 significance threshold
    sig_threshold = -np.log10(sig_threshold)
    sig_var = f'-log10({sig_var})'

    # Plot volcano plot
    fig = volcano_plot(
        data=data_to_plot,
        effect_size_var=enrichment_var,
        sig_var=sig_var,
        sig_threshold=sig_threshold,
        high_effect_size_threshold=enrichment_threshold,
        hover_data=hover_data,
    )

    return fig


# Callback for motif enrichment table
@callback(
    Output('enriched-motifs-table', 'children'),
    [
        Input('run-selector', 'value'), 
        Input('program-selector', 'value'), 
        Input('enriched-motifs-threshold', 'value'),
        Input('enriched-motifs-enrichment-threshold', 'value'),
        Input('table-motifs-E_P_type-selector', 'value'),
        Input('table-motifs-database-selector', 'value'),
        Input('table-motifs-test_type-selector', 'value'),
        Input('table-motifs-level_key-selector', 'value')
    ]
)
def update_enriched_motifs_table(
    selected_run, 
    selected_program, 
    sig_threshold,
    enrichment_threshold,
    e_p_type,
    database,
    test_type,
    level_key,
    categorical_var = "program_name",
    sig_var = "adj_pval",
    enrichment_var = "stat",
    debug=False,
):
  
    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program}, sig_threshold: {sig_threshold}")
        print(f"E/P type: {e_p_type}, Database: {database}, Test type: {test_type}, Level key: {level_key}")
        print(f"Enrichment threshold: {enrichment_threshold}")

    # Retrieve the relevant data
    data_to_plot = results['motif_enrichments'][selected_run]["results"][f"{e_p_type}_{database}_{test_type}_{motif_enrichment_stratification_key}_{level_key}"]

    # Make sure x-axis is string
    data_to_plot[categorical_var] = data_to_plot[categorical_var].astype(str)

    # Filter data for the selected program
    data_to_plot = data_to_plot[data_to_plot[categorical_var] == selected_program]

    # Filter data for significant terms
    data_to_plot = data_to_plot[(data_to_plot[sig_var] < sig_threshold) & (data_to_plot[enrichment_var] > enrichment_threshold)]
    
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
    [
        Input('run-selector', 'value'), 
        Input('program-selector', 'value'),
        Input('table-motifs-E_P_type-selector', 'value'),
        Input('table-motifs-database-selector', 'value'),
        Input('table-motifs-test_type-selector', 'value'),
        Input('table-motifs-level_key-selector', 'value'),
        Input('enriched-motifs-threshold', 'value'),
        Input('enriched-motifs-enrichment-threshold', 'value'),
    ]
)
def update_volcano_plot_motifs(
    selected_run,
    selected_program,
    e_p_type,
    database,
    test_type,
    level_key,
    sig_threshold,
    enrichment_threshold,
    categorical_var = "motif",
    sig_var = "adj_pval",
    enrichment_var = "stat",
    hover_data=["program_name", "motif", "adj_pval", "stat"],
    debug=False,
):

    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program}")
        print(f"E/P type: {e_p_type}, Database: {database}, Test type: {test_type}, Level key: {level_key}")
        print(f"Sig threshold: {sig_threshold}, Enrichment threshold: {enrichment_threshold}")

    # Assuming we want to plot something from the selected run
    data_to_plot = results['motif_enrichments'][selected_run]["results"][f"{e_p_type}_{database}_{test_type}_{motif_enrichment_stratification_key}_{level_key}"]

    # Filter data for the selected program
    data_to_plot = data_to_plot.query(f"program_name == '{selected_program}'")

    # -log10 transform adj_pvals
    data_to_plot[f'-log10({sig_var})'] = -np.log10(data_to_plot[sig_var])

    # Make sure x-axis is string
    data_to_plot[categorical_var] = data_to_plot[categorical_var].astype(str)

    # -log10 significance threshold
    sig_treshold = -np.log10(sig_threshold)
    sig_var = f'-log10({sig_var})'

    # Plot volcano plot
    fig = volcano_plot(
        data=data_to_plot,
        effect_size_var=enrichment_var,
        sig_var=sig_var,
        sig_threshold=sig_treshold,
        high_effect_size_threshold=enrichment_threshold,
        hover_data=hover_data,
    )

    return fig


# Callback for trait enrichment table
@callback(
    Output('enriched-traits-table', 'children'),
    [
        Input('run-selector', 'value'), 
        Input('program-selector', 'value'), 
        Input('enriched-traits-threshold', 'value'),
        Input('enriched-traits-enrichment-threshold', 'value'),
        Input('table-traits-library-selector', 'value'),
        Input('table-traits-method-selector', 'value'),
        Input('trait-type-selector', 'value'),
        Input('trait-category-selector', 'value')
    ]
)
def update_enriched_traits_table(
    selected_run, 
    selected_program,
    sig_threshold,
    enrichment_threshold,
    library,
    method,
    trait_type,
    trait_categories,
    categorical_var = "program_name",
    sig_var = "adj_pval",
    enrichment_var = "enrichment",
    debug=False,
):

    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program}, sig_threshold: {sig_threshold}, enrichment_threshold: {enrichment_threshold}")
        print(f"Library: {library}, Method: {method}, Trait Type: {trait_type}")

    # Retrieve the relevant data
    data_to_plot = results['trait_enrichments'][selected_run]["results"][f"{library}_{method}"]

    # Make sure x-axis is string
    data_to_plot[categorical_var] = data_to_plot[categorical_var].astype(str)

    # Filter data for the selected program
    data_to_plot = data_to_plot[data_to_plot[categorical_var] == selected_program]

    # Filter data for significant terms
    data_to_plot = data_to_plot[(data_to_plot[sig_var] < sig_threshold) & (data_to_plot[enrichment_var] > enrichment_threshold)]

    # Filter data for binary or continuous traits
    if trait_type == 'binary':
        data_to_plot = data_to_plot.query("trait_category != 'measurement'")
    elif trait_type == 'continuous':
        data_to_plot = data_to_plot.query("trait_category == 'measurement'")
    else:
        print("Invalid trait type.")
        return html.Div("Invalid trait type.")
    
    # Filter data for selected trait categories
    data_to_plot = data_to_plot[data_to_plot["trait_category"].isin(trait_categories)]

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


# Callback for phewas binary plot
@callback(
    Output('phewas-program-plot', 'figure'),
    [
        Input('run-selector', 'value'), 
        Input('program-selector', 'value'),
        Input('table-traits-library-selector', 'value'),
        Input('table-traits-method-selector', 'value'),
        Input('trait-type-selector', 'value'),
        Input('trait-category-selector', 'value'),
        Input('enriched-traits-threshold', 'value'),
        Input('enriched-traits-enrichment-threshold', 'value'),
    ]
)
def update_phewas_plot(
    selected_run,
    selected_program,
    library,
    method,
    trait_type,
    trait_categories,
    sig_threshold,
    enrichment_threshold,
    debug=False,
):
    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program}, trait_type: {trait_type}")
        print(f"Library: {library}, Method: {method}")

    # Assuming we want to plot something from the selected run
    data_to_plot = results['trait_enrichments'][selected_run]["results"][f"{library}_{method}"]

    # Make sure x-axis is string
    data_to_plot["program_name"] = data_to_plot["program_name"].astype(str)

    # Filter data for the selected program
    data_to_plot = data_to_plot.query(f"program_name == '{selected_program}'")

    # Filter data for the selected trait type
    if trait_type == 'binary':
        data_to_plot = data_to_plot.query("trait_category != 'measurement'")
    elif trait_type == 'continuous':
        data_to_plot = data_to_plot.query("trait_category == 'measurement'")
    else:
        print("Invalid trait type.")
        return go.Figure()
    
    # Filter data for selected trait categories
    data_to_plot = data_to_plot[data_to_plot["trait_category"].isin(trait_categories)]
    
    # Plot
    fig = px.scatter(
        data_to_plot,
        x='trait_reported',
        y='-log10(adj_pval)',
        color='trait_category',
        title="",
        hover_data=["program_name", "trait_reported", "trait_category", "adj_pval", "study_id", "pmid"]
    )

    # Customize layout
    fig.update_layout(
        xaxis_title='trait_reported',
        yaxis_title='-log10(adj_pval)',
        yaxis=dict(tickformat=".1f"),
        width=1000,
        height=800,
        xaxis=dict(showticklabels=False),
        plot_bgcolor='rgba(0,0,0,0)',  # Transparent background
    )

    # Add horizontal dashed line for significance threshold
    fig.add_hline(
        y=-np.log10(sig_threshold),
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
    
    curr_obs = results['obs']
    curr_obs_membership = results['obs_memberships'][selected_run]
    curr_obs_membership.index = curr_obs.index
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
    Output('program-dim-reduction-scatter', 'src'),
    [Input('program-dim-reduction-selector', 'value'), Input('run-selector', 'value'), Input('program-selector', 'value')]
)
def update_program_dim_reduction_plot(
    selected_dim_reduction,
    selected_run,
    selected_program,
    size=1,
    debug=False,
    static=True
):
    
    if debug:
        print(f"Selected dim reduction: {selected_dim_reduction}, selected run: {selected_run}, selected program: {selected_program}")
        print(f"Static: {static}")

    if selected_dim_reduction is None:
        return go.Figure(), html.Div()  # Return an empty figure and legend placeholder

    # Get vector of continuous values for the program membership
    curr_obs = results['obs']
    curr_obs_membership = results['obs_memberships'][selected_run][selected_program]

    # Extract the relevant data for plotting (they will be selected_dim_reduction_0 and selected_dim_reduction_1)
    data = obsms[selected_dim_reduction][[f'{selected_dim_reduction}_0', f'{selected_dim_reduction}_1']]
    data.columns = ['X', 'Y']

    colors = curr_obs_membership
    categorical_color_map = None
    continuous_color_map = "viridis"
    if debug:
        print(f"Colors: {colors[:5]}")
        print(f"Color map: {continuous_color_map}")

    if static:

        if debug:
            print("Generating static plot using matplotlib")

        # Plot using static matplotlib function
        fig, ax = plt.subplots()
        scatterplot_static(
            ax=ax,
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

        cbar = plt.colorbar(ax.collections[0], ax=ax)
        cbar.set_label("Program Loadings")
    
        # Save the figure to a temporary file
        fig = fig_to_uri(fig)

        return fig

    else:
        raise NotImplementedError("Interactive scatter plot not implemented yet.")


# Callback for perturbation table
@callback(
    Output('perturbation-table', 'children'),
    [
        Input('run-selector', 'value'), 
        Input('program-selector', 'value'), 
        Input('perturbation-threshold', 'value'),
        Input('table-perturbations-gene_guide-selector', 'value'),
        Input('table-perturbations-level_key-selector', 'value')
    ]
)
def update_perturbation_table(
    selected_run, 
    selected_program, 
    sig_threshold,
    gene_guide,
    level_key,
    categorical_var = "target_name",
    sig_var = "adj_pval",
    debug=False,
):
    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program}, sig_threshold: {sig_threshold}")
        print(f"Gene guide: {gene_guide}, Level key: {level_key}")

    # Retrieve the relevant data
    data_to_plot = results['perturbation_associations'][selected_run]["results"][f"{gene_guide}_{perturbation_association_stratification_key}_{level_key}"]

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
    [
        Input('run-selector', 'value'), 
        Input('program-selector', 'value'),
        Input('table-perturbations-gene_guide-selector', 'value'),
        Input('table-perturbations-level_key-selector', 'value'),
    ]
)
def update_perturbation_association_plot(
    selected_run, 
    selected_program,
    selected_gene_guide,
    selected_level_key,
    categorical_var = "target_name",
    sig_var = "adj_pval",
    sig_threshold = 0.05,
    effect_size_var = "log2FC",
    low_effect_size_threshold = -0.5,
    high_effect_size_threshold = 0.5,
    debug=False
):
    
    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program} for perturbation association plot")
        print(f"Gene guide: {selected_gene_guide}, Level key: {selected_level_key}")

    # Assuming we want to plot something from the selected run
    data_to_plot = results['perturbation_associations'][selected_run]["results"][f"{selected_gene_guide}_{perturbation_association_stratification_key}_{selected_level_key}"]

    # Make sure x-axis is string
    data_to_plot["program_name"] = data_to_plot["program_name"].astype(str)

    # Filter data for the selected program
    data_to_plot = data_to_plot.query(f"program_name == '{selected_program}'")

    # -log10 transform adj_pvals
    data_to_plot[f'-log10({sig_var})'] = -np.log10(data_to_plot[sig_var])

    # Make sure x-axis is string
    data_to_plot[categorical_var] = data_to_plot[categorical_var].astype(str)

    # -log10 significance threshold
    sig_threshold = -np.log10(sig_threshold)

    print(f"data_to_plot: {data_to_plot.head()}")

    # Plot volcano plot
    fig = volcano_plot(
        data=data_to_plot,
        effect_size_var=effect_size_var,
        sig_var=f'-log10({sig_var})',
        sig_threshold=sig_threshold,
        low_effect_size_threshold=low_effect_size_threshold,
        high_effect_size_threshold=high_effect_size_threshold,
        hover_data=["program_name", "target_name", "adj_pval", "stat"],
    )

    return fig


# Callback to save annotations
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
    annotations_file = os.path.join(path_report_out, annotations_loc)
    if os.path.exists(annotations_file):
        df = pd.read_csv(annotations_file)
    else:
        df = pd.DataFrame(columns=['Program', 'Annotation', 'TimeStamp'])
    
    # Append the annotation
    new_row = {'Program': selected_program, 'Annotation': annotation_text, 'TimeStamp': datetime.datetime.now()}
    df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
    
    # Save the updated DataFrame back to the CSV file
    df.to_csv(annotations_file, index=False)
    
    return f"Annotation for '{selected_program}' has been saved."
