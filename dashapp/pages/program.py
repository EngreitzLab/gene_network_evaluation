import os
import dash
import pickle
import pandas as pd
import dash_bootstrap_components as dbc
from dash import html, dcc, callback, Input, Output
from dash import dash_table
from plot import scatterplot, barplot, lollipop_plot, boxplot

# Ouput directory
path_pipeline_outs = "/cellar/users/aklie/opt/gene_program_evaluation/dashapp/example_data/iPSC_EC_evaluations"
data_key = "rna"

dash.register_page(__name__, order=3)

# Load results from pickle
with open(os.path.join(path_pipeline_outs, "results.pkl"), "rb") as f:
    results = pickle.load(f)

# Get the first data type and its corresponding second-level keys as default
default_data_type = "explained_variance_ratios"
default_run = list(results[default_data_type].keys())[0]

# Get the first program from the data type and run as default
programs = sorted(list(results[default_data_type][default_run]["program_name"].astype(str).unique()))
default_program = programs[0]

layout = dbc.Container([
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
            dcc.Graph(id='gene-loadings-plot'),
        ]),

        # Tab 2: Gene Set Enrichment
        dcc.Tab(label='Gene Set Enrichment', children=[
            html.H2("Gene Set Enrichment Results", className="mt-4 mb-3"),
            
            dbc.Row([
                dbc.Col([
                    html.H3("Table of Enriched Terms"),
                    html.Div(id='enriched-terms-table', className="mb-4"),
                ], width=12),
            ]),

            html.H3("GSEA Plot for Selected Term", className="mt-3 mb-3"),
            dcc.Graph(id='gsea-plot'),
        ]),

        # Tab 3: Motif Enrichment
        dcc.Tab(label='Motif Enrichment', children=[
            html.H2("Motif Enrichment Results", className="mt-4 mb-3"),
            dbc.Row([
                dbc.Col([
                    html.H3("Table of Enriched Motifs"),
                    html.Div(id='enriched-motifs-table', className="mb-4"),
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
                    html.Div(id='enriched-traits-table', className="mb-4"),
                ], width=12),
            ]),
            
            html.H3("Phewas Plot for Selected Program", className="mt-3 mb-3"),
            dcc.Graph(id='phewas-plot'),
        ]),

        # Tab 5: Categorical Association
        dcc.Tab(label='Categorical Association', children=[
            html.H2("Categorical Association Analysis", className="mt-4 mb-3"),
            html.P("Boxplot", className="mb-3"),
            dcc.Graph(id='covariate-association-plot'),
        ]),

        # Tab 6: Perturbation
        dcc.Tab(label='Perturbation', children=[
            html.H2("Perturbation", className="mt-4 mb-3"),
            html.P("Add any perturbation analysis here.", className="mb-3"),
            dcc.Graph(id='perturbation-plot'),
        ]),

        # Tab 7: Annotations
        dcc.Tab(label='Annotations', children=[
            html.H2("Type in Annotation", className="mt-4 mb-3"),
            dbc.Row([
                dbc.Col([
                    dcc.Textarea(
                        id='annotation-text',
                        placeholder="Enter annotations or notes here...",
                        style={'width': '100%', 'height': 100},
                    )
                ], width=12),
            ], className="mb-4"),
        ]),
    ])
], fluid=True, className="p-4")


@callback(
    Output('gene-loadings-plot', 'figure'),
    [Input('run-selector', 'value'), Input('program-selector', 'value')]
)
def update_plot(selected_run, selected_program, debug=True):
    n = 25
    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program}")

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
    [Input('run-selector', 'value'), Input('program-selector', 'value')]
)
def update_enriched_terms_table(selected_run, selected_program, debug=True):
    categorical_var = "program_name"
    count_var = "Term"
    sig_var = "FDR q-val"
    sig_threshold = 0.25

    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program}")

    # Retrieve the relevant data
    data_to_plot = results['geneset_enrichments'].get(selected_run, pd.DataFrame())

    # Make sure x-axis is string
    data_to_plot['program_name'] = data_to_plot['program_name'].astype(str)

    if data_to_plot.empty:
        print("No data found for the selected run.")
        return html.Div("No data available for the selected run and program.")

    # Filter data for the selected program
    data_to_plot = data_to_plot[data_to_plot[categorical_var] == selected_program]

    if data_to_plot.empty:
        print("No data found for the selected program.")
        return html.Div("No enriched terms available for this program.")

    # Filter data for significant terms
    data_to_plot = data_to_plot[data_to_plot[sig_var] < sig_threshold]

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
        style_data_conditional=[
            {
                'if': {'filter_query': f'{{{sig_var}}} < {sig_threshold}', 'column_id': sig_var},
                'backgroundColor': 'rgb(230, 240, 250)',
                'color': 'black',
            },
        ],
    )

    return table

# Callback for motif enrichment table
@callback(
    Output('enriched-motifs-table', 'children'),
    [Input('run-selector', 'value'), Input('program-selector', 'value')]
)
def update_enriched_motifs_table(selected_run, selected_program, debug=True):
    categorical_var = "program_name"
    count_var = "motif"
    sig_var = "pval"
    sig_threshold = 0.75

    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program}")

    # Retrieve the relevant data
    data_to_plot = results['motif_enrichments'].get(selected_run, pd.DataFrame())

    # Make sure x-axis is string
    data_to_plot['program_name'] = data_to_plot['program_name'].astype(str)

    if data_to_plot.empty:
        print("No data found for the selected run.")
        return html.Div("No data available for the selected run and program.")

    # Filter data for the selected program
    data_to_plot = data_to_plot[data_to_plot[categorical_var] == selected_program]

    if data_to_plot.empty:
        print("No data found for the selected program.")
        return html.Div("No enriched motifs available for this program.")

    # Filter data for significant terms
    data_to_plot = data_to_plot[data_to_plot[sig_var] < sig_threshold]

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
        style_data_conditional=[
            {
                'if': {'filter_query': f'{{{sig_var}}} < {sig_threshold}', 'column_id': sig_var},
                'backgroundColor': 'rgb(230, 240, 250)',
                'color': 'black',
            },
        ],
    )

    return table


# Callback for trait enrichment table
@callback(
    Output('enriched-traits-table', 'children'),
    [Input('run-selector', 'value'), Input('program-selector', 'value')]
)
def update_enriched_traits_table(selected_run, selected_program, debug=True):
    categorical_var = "program_name"
    count_var = "Term"
    sig_var = "FDR q-val"
    sig_threshold = 0.25

    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program}")

    # Retrieve the relevant data
    data_to_plot = results['trait_enrichments'].get(selected_run, pd.DataFrame())

    # Make sure x-axis is string
    data_to_plot['program_name'] = data_to_plot['program_name'].astype(str)

    if data_to_plot.empty:
        print("No data found for the selected run.")
        return html.Div("No data available for the selected run and program.")

    # Filter data for the selected program
    data_to_plot = data_to_plot[data_to_plot[categorical_var] == selected_program]

    if data_to_plot.empty:
        print("No data found for the selected program.")
        return html.Div("No enriched traits available for this program.")

    # Filter data for significant terms
    data_to_plot = data_to_plot[data_to_plot[sig_var] < sig_threshold]

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
        style_data_conditional=[
            {
                'if': {'filter_query': f'{{{sig_var}}} < {sig_threshold}', 'column_id': sig_var},
                'backgroundColor': 'rgb(230, 240, 250)',
                'color': 'black',
            },
        ],
    )

    return table


# Callback for covariate association plot
@callback(
    Output('covariate-association-plot', 'figure'),
    [Input('run-selector', 'value'), Input('program-selector', 'value')]
)
def update_covariate_association_plot(selected_run, selected_program, debug=True):

    if debug:
        print(f"Selected run: {selected_run}, program: {selected_program}")

    curr_obs = results['obs'][selected_run]
    curr_obs_membership = results['obs_memberships'][selected_run]
    concat_df = pd.concat([curr_obs_membership[selected_program], curr_obs["sample"]], axis=1)
    
    fig = boxplot(
        concat_df, 
        x_column=selected_program, 
        y_column="sample", 
        title="",
        x_axis_title=selected_program,
        y_axis_title="Sample",
    )

    return fig