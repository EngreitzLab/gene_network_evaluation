import os
import dash
import pickle
from dash import html, dcc, callback, Input, Output
import dash_bootstrap_components as dbc
import plotly.express as px
import pandas as pd
from plot import scatterplot, barplot
from utils import count, count_unique

# Ouput directory
path_pipeline_outs = "/cellar/users/aklie/opt/gene_program_evaluation/dashapp/example_data/iPSC_EC_evaluations"
phewas_metadata = "/cellar/users/aklie/opt/gene_program_evaluation/smk/resources/OpenTargets_L2G_Filtered.csv.gz"
data_key = "rna"

dash.register_page(__name__, order=2)

# Load results from pickle
with open(os.path.join(path_pipeline_outs, "results.pkl"), "rb") as f:
    results = pickle.load(f)

# Get the first data type and its corresponding second-level keys as default
default_data_type = list(results.keys())[0]
default_run = list(results[default_data_type].keys())[0]

layout = dbc.Container([
    html.H1("Single Run Analysis", className="mb-4"),

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

    # Tabs for different sections
    dcc.Tabs([
        dcc.Tab(label='Goodness of Fit', children=[
            html.H2("Goodness of Fit", className="mt-3 mb-3"),
            
            dbc.Row([
                dbc.Col([
                    html.H3("Explained Variance Per Component"),
                    dcc.Graph(id='explained-variance-plot'),
                ], width=6),

                dbc.Col([
                    html.H3("Unique Genesets Enriched"),
                    dcc.Graph(id='num-enriched-genesets'),
                ], width=6),
            ], className="mb-4"),

            dbc.Row([
                dbc.Col([
                    html.H3("Unique Motifs Enriched"),
                    dcc.Graph(id='num-enriched-motifs'),
                ], width=6),

                dbc.Col([
                    html.H3("Unique Traits Enriched"),
                    dcc.Graph(id='num-enriched-traits'),
                ], width=6),
            ], className="mb-4"),

            html.H3("Perturbation Associations", className="mt-4"),
            dcc.Graph(id='perturbation-associations'),
        ]),

        dcc.Tab(label='Covariate Association', children=[
            html.H2("Covariate Association", className="mt-3 mb-3"),

            dbc.Row([
                dbc.Col([
                    html.H3("Topic-Trait Correlation Heatmap"),
                    dcc.Graph(id='topic-trait-heatmap'),
                ], width=12),
            ], className="mb-4"),

            dbc.Row([
                dbc.Col([
                    html.H3("Structure Plot"),
                    dcc.Graph(id='structure-plot'),
                ], width=12),
            ], className="mb-4"),

            dbc.Row([
                dbc.Col([
                    html.H3("Enrichment with Respect to Covariate X"),
                    dcc.Graph(id='covariate-enrichment'),
                ], width=12),
            ], className="mb-4"),
        ]),

        dcc.Tab(label='Perturbation Association', children=[
            html.H2("Perturbation Association", className="mt-3 mb-3"),
        ]),

        dcc.Tab(label='Trait Enrichment', children=[
            html.H2("Trait Enrichment", className="mt-3 mb-3"),

            dbc.Row([
                dbc.Col([
                    html.H3("Phewas Plot"),
                    dcc.Graph(id='phewas-plot-binary'),
                ], width=12),
            ], className="mb-4"),
        ]),
    ])
], fluid=True, className="p-4")


@callback(
    Output('explained-variance-plot', 'figure'),
    [Input('run-selector', 'value')]
)
def update_plot(selected_run, debug=True):

    if debug:
        print(f"Selected run: {selected_run}")

    # Assuming we want to plot something from the selected run
    data_to_plot = results['explained_variance_ratios'][selected_run]  # Adjust to your actual plot logic

    # Make sure x-axis is string
    data_to_plot['program_name'] = data_to_plot['program_name'].astype(str)

    # Example plot using Plotly (replace with your actual plotting code)
    fig = scatterplot(
        data=data_to_plot,
        x_column='program_name',
        y_column='explained_variance_ratio',
        title='',
        x_axis_title='Component',
        y_axis_title='Variance Explained (R²)',
    )
    return fig


@callback(
    Output('num-enriched-genesets', 'figure'),
    [Input('run-selector', 'value')]
)
def update_plot(selected_run, debug=True):
    categorical_var = "program_name"
    count_var = "Term"
    sig_var = "FDR q-val"
    sig_threshold = 0.25

    if debug:
        print(f"Selected run: {selected_run}")

    # Assuming we want to plot something from the selected run
    data_to_plot = results['geneset_enrichments'][selected_run]  # Adjust to your actual plot logic

    # Make sure x-axis is string
    data_to_plot['program_name'] = data_to_plot['program_name'].astype(str)

    filtered_data = data_to_plot[data_to_plot[sig_var] < sig_threshold]

    # Get the count of all gene sets passing the FDR q-value cutoff for each program
    count_df = count(categorical_var=categorical_var, count_var=count_var, dataframe=filtered_data)

    # Get the count of unique gene sets passing the FDR q-value cutoff for each program
    unique_data = filtered_data.sort_values(by=sig_var)
    unique_data = unique_data.drop_duplicates(subset=count_var)
    unique_df = count_unique(categorical_var=categorical_var, count_var=count_var, dataframe=unique_data)
    unique_df = unique_df.sort_values(count_var, ascending=False)

    # Example plot using Plotly (replace with your actual plotting code)
    fig = barplot(
        data=unique_df,
        x_column='program_name',
        y_column='Term',
        title='',
        x_axis_title='Program',
        y_axis_title='Count',
        show_xaxis_labels=False
    )
    return fig

@callback(
    Output('num-enriched-motifs', 'figure'),
    [Input('run-selector', 'value')]
)
def update_plot(selected_run, debug=True):
    categorical_var = "program_name"
    count_var = "motif"
    sig_var = "pval"
    sig_threshold = 0.75

    if debug:
        print(f"Selected run: {selected_run}")

    # Assuming we want to plot something from the selected run
    data_to_plot = results['motif_enrichments'][selected_run]  # Adjust to your actual plot logic

    # Make sure x-axis is string
    data_to_plot['program_name'] = data_to_plot['program_name'].astype(str)

    filtered_data = data_to_plot[data_to_plot[sig_var] < sig_threshold]

    # Get the count of all gene sets passing the FDR q-value cutoff for each program
    count_df = count(categorical_var=categorical_var, count_var=count_var, dataframe=filtered_data)

    # Get the count of unique gene sets passing the FDR q-value cutoff for each program
    unique_data = filtered_data.sort_values(by=sig_var)
    unique_data = unique_data.drop_duplicates(subset=count_var)
    unique_df = count_unique(categorical_var=categorical_var, count_var=count_var, dataframe=unique_data)
    unique_df = unique_df.sort_values(count_var, ascending=False)

    # Example plot using Plotly (replace with your actual plotting code)
    fig = barplot(
        data=count_df,
        x_column=categorical_var,
        y_column=count_var,
        title='',
        x_axis_title='Program',
        y_axis_title='Count',
        show_xaxis_labels=False
    )
    return fig


@callback(
    Output('num-enriched-traits', 'figure'),
    [Input('run-selector', 'value')]
)
def update_plot(selected_run, debug=True):
    categorical_var = "program_name"
    count_var = "Term"
    sig_var = "FDR q-val"
    sig_threshold = 0.25

    if debug:
        print(f"Selected run: {selected_run}")

    # Assuming we want to plot something from the selected run
    data_to_plot = results['trait_enrichments'][selected_run]

    # Make sure x-axis is string
    data_to_plot['program_name'] = data_to_plot['program_name'].astype(str)

    filtered_data = data_to_plot[data_to_plot[sig_var] < sig_threshold]
    
    # Get the count of all gene sets passing the FDR q-value cutoff for each program
    count_df = count(categorical_var=categorical_var, count_var=count_var, dataframe=filtered_data)

    # Get the count of unique gene sets passing the FDR q-value cutoff for each program
    unique_data = filtered_data.sort_values(by=sig_var)
    unique_data = unique_data.drop_duplicates(subset=count_var)
    unique_df = count_unique(categorical_var=categorical_var, count_var=count_var, dataframe=unique_data)
    unique_df = unique_df.sort_values(count_var, ascending=False)

    # Example plot using Plotly (replace with your actual plotting code)
    fig = barplot(
        data=count_df,
        x_column=categorical_var,
        y_column=count_var,
        title='',
        x_axis_title='Program',
        y_axis_title='Count',
        show_xaxis_labels=False
    )
    return fig


@callback(
    Output('perturbation-associations', 'figure'),
    [Input('run-selector', 'value')]
)
def update_plot(selected_run, debug=True):
    categorical_var = "program"
    count_var = "guide_name"
    sig_var = "pval"
    sig_threshold = 0.25

    if debug:
        print(f"Selected run: {selected_run}")

    # Assuming we want to plot something from the selected run
    data_to_plot = results['perturbation_associations'][selected_run]

    # Make sure x-axis is string
    data_to_plot[categorical_var] = data_to_plot[categorical_var].astype(str)

    filtered_data = data_to_plot[data_to_plot[sig_var] < sig_threshold]

    # Get the count of all gene sets passing the FDR q-value cutoff for each program
    count_df = count(categorical_var=categorical_var, count_var=count_var, dataframe=filtered_data)

    # Get the count of unique gene sets passing the FDR q-value cutoff for each program
    unique_data = filtered_data.sort_values(by=sig_var)
    unique_data = unique_data.drop_duplicates(subset=count_var)
    unique_df = count_unique(categorical_var=categorical_var, count_var=count_var, dataframe=unique_data)
    unique_df = unique_df.sort_values(count_var, ascending=False)

    # Example plot using Plotly (replace with your actual plotting code)
    fig = barplot(
        data=count_df,
        x_column=categorical_var,
        y_column=count_var,
        title='',
        x_axis_title='Program',
        y_axis_title='Count',
        show_xaxis_labels=False
    )
    return fig


import numpy as np
def process_enrichment_data(enrich_res,
                            metadata,
                            pval_col="FDR q-val",
                            enrich_geneset_id_col="Term",
                            metadata_geneset_id_col="trait_efos",
                            color_category_col="trait_category",
                            program_name_col="program_name",
                            annotation_cols=["trait_reported", "Lead_genes", "study_id", "pmid"]):

    # Read in enrichment results
    if isinstance(enrich_res, str):
        enrich_df = pd.read_csv(enrich_res)
    elif isinstance(enrich_res, pd.DataFrame):
        enrich_df = enrich_res
    else:
        raise ValueError("enrich_res must be either a pandas DataFrame or a file path to a CSV file.")

    if isinstance(metadata, str):
        metadata_df = pd.read_csv(metadata, compression='gzip', low_memory=False)
    elif isinstance(metadata, pd.DataFrame):
        metadata_df = metadata
    else:
        raise ValueError("metadata must be either a pandas DataFrame or a file path to a CSV file.")

    # Join the enrichment results and the metadata
    enrich_ps = enrich_df.merge(metadata_df, left_on=enrich_geneset_id_col, right_on=metadata_geneset_id_col, how="left")

    # Only keep the relevant columns
    keep_cols = list([enrich_geneset_id_col, pval_col, metadata_geneset_id_col, color_category_col, program_name_col] + annotation_cols)
    enrich_ps = enrich_ps[keep_cols]

    # Sort by P-value
    enrich_ps = enrich_ps.drop_duplicates().sort_values(by=[color_category_col, pval_col])

    # If the input P-value == 0, then replace it with the lowest non-zero P-value in the dataframe
    min_value = enrich_ps.query(f"`{pval_col}` > 0")[pval_col].min()

    # Compute the -log(10) P-value and deal with edge-cases (e.g. P=0, P=1)
    enrich_ps.loc[enrich_ps[pval_col] == 0, pval_col] = min_value  # Replace P=0 with min non-0 p-value
    enrich_ps['-log10(p-value)'] = abs(-1 * np.log10(enrich_ps[pval_col]))

    enrich_ps.reset_index(drop=True, inplace=True)

    return enrich_ps


@callback(
    Output('phewas-plot-binary', 'figure'),
    Input('run-selector', 'value')
)
def update_phewas_plot_binary(selected_run):

    # Assuming we want to plot something from the selected run
    data_to_plot = results['trait_enrichments'][selected_run]
    data_to_plot = process_enrichment_data(
        enrich_res=data_to_plot,
        metadata=phewas_metadata,
    )

    fig = px.scatter(
        data_to_plot.query("trait_category != 'measurement'"),
        x='trait_reported',
        y='-log10(p-value)',
        color='trait_category',
        title="Endothelial Cell Programs x GWAS Binary Outcome Enrichments",
        hover_data=["program_name", "trait_reported", "trait_category", "FDR q-val", "Lead_genes", "study_id", "pmid"]
    )

    # Customize layout
    fig.update_layout(
        xaxis_title='trait_reported',
        yaxis_title='-log10(p-value)',
        yaxis=dict(tickformat=".1f"),
        width=1000,
        height=800,
        xaxis_tickfont=dict(size=4)
    )

    # Add horizontal dashed line for significance threshold
    fig.add_hline(y=-np.log10(0.05), line_dash="dash",
                    annotation_text='Significance Threshold (0.05)', annotation_position="top right")

    return fig
