import os
import dash
import pickle
from dash import html, dcc, callback, Input, Output
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import pandas as pd
from plot import scatterplot, barplot
from utils import filter_and_count
import diskcache


# Register the page
dash.register_page(__name__, order=1)

# Get the app and cache
app = dash.get_app()
cache = diskcache.Cache("./.cache")
results = cache.get("results")

# Get the first data type and its corresponding second-level keys as default
default_data_type = list(results.keys())[0]
default_run = list(results[default_data_type].keys())[0]

# Create layout
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

        # 1. Tab for Goodness of Fit
        dcc.Tab(label='Goodness of Fit', children=[
            html.H2("Goodness of Fit", className="mt-3 mb-3"),
            
            # Explained Variance Per Program and Number of Genesets Enriched Per Program
            dbc.Row([
                dbc.Col([
                    html.H3("Explained Variance Per Program"),
                    dcc.Graph(id='explained-variance-plot'),
                ], width=6),

                dbc.Col([
                    html.H3("Number of Genesets Enriched Per Program"),
                    dcc.Graph(id='num-enriched-genesets'),
                ], width=6),
            ], className="mb-4"),

            # Number of Motifs Enriched Per Program and Number of Traits Enriched Per Program
            dbc.Row([
                dbc.Col([
                    html.H3("Number of Motifs Enriched Per Program"),
                    dcc.Graph(id='num-enriched-motifs'),
                ], width=6),

                dbc.Col([
                    html.H3("Number of Traits Enriched Per Program"),
                    dcc.Graph(id='num-enriched-traits'),
                ], width=6),
            ], className="mb-4"),

            # Number of Perturbation Associations Per Program
            dbc.Row([
                dbc.Col([
                    html.H3("Number of Perturbation Associations Per Program", className="mt-4"),
                    dcc.Graph(id='num-perturbation-associations'),
                ], width=6),
            ], className="mb-4"),
        ]),

        # 2. Tab for Covariate Association
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

        dcc.Tab(label='Trait Enrichment', children=[
            html.H2("Trait Enrichment", className="mt-3 mb-3"),

            dbc.Row([
                dbc.Col([
                    html.H3("GEP x GWAS Binary Outcome Enrichments"),
                    dcc.Graph(id='phewas-binary-plot'),
                ], width=12),
            ], className="mb-4"),

            dbc.Row([
                dbc.Col([
                    html.H3("GEP x GWAS Continuous Outcome Enrichments"),
                    dcc.Graph(id='phewas-continuous-plot'),
                ], width=12),
            ], className="mb-4"),
        ]),
    ])
], fluid=True, className="p-4")


@callback(
    Output('explained-variance-plot', 'figure'),
    [Input('run-selector', 'value')]
)
def update_explained_variance_plot(
    selected_run, 
    debug=False,
    categorical_var = "program_name",
    ):

    if debug:
        print(f"Selected run for explained variance plot: {selected_run}")

    # Assuming we want to plot something from the selected run
    data_to_plot = results['explained_variance_ratios'][selected_run]  # Adjust to your actual plot logic

    # Make sure x-axis is string
    data_to_plot[categorical_var] = data_to_plot[categorical_var].astype(str)

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
def update_num_enriched_genesets_plot(
    selected_run, 
    debug=True,
    unique=False,
    categorical_var = "program_name",
    count_var = "term",
    sig_var = "adj_pval",
    sig_threshold = 0.05
):

    if debug:
        print(f"Selected run for gene set enrichment plot: {selected_run}")
        print(f"Filtering gene set enrichments with {sig_var} < {sig_threshold}")
        print(f"Counting unique {count_var} enriched for each {categorical_var}")

    # Assuming we want to plot something from the selected run
    data_to_plot = results['geneset_enrichments'][selected_run]  # Adjust to your actual plot logic

    # Make sure x-axis is string
    data_to_plot[categorical_var] = data_to_plot[categorical_var].astype(str)

    # Filter data
    count_df, unique_df = filter_and_count(
        data=data_to_plot,
        categorical_var=categorical_var,
        count_var=count_var,
        sig_var=sig_var,
        sig_threshold=sig_threshold
    )

    # Choose whether to plot all gene sets or only unique gene sets
    if unique:
        data_to_plot = unique_df
    else:
        data_to_plot = count_df

    if debug:
        print(f"Counts of gene sets enriched: {data_to_plot}")

    # Example plot using Plotly (replace with your actual plotting code)
    fig = go.Figure(layout=dict(template='plotly'))
    fig = barplot(
        data=data_to_plot,
        x_column=categorical_var,
        y_column=count_var,
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
def update_num_enriched_motifs_plot(
    selected_run, 
    debug=True,
    unique=False,
    categorical_var = "program_name",
    count_var = "motif",
    sig_var = "pval",
    sig_threshold = 0.05
):

    if debug:
        print(f"Selected run for motif enrichment plot: {selected_run}")
        print(f"Filtering motif enrichments with {sig_var} < {sig_threshold}")
        print(f"Counting unique {count_var} enriched for each {categorical_var}")
    
    # Assuming we want to plot something from the selected run
    data_to_plot = results['motif_enrichments'][selected_run]  # Adjust to your actual plot logic

    # Make sure x-axis is string
    data_to_plot[categorical_var] = data_to_plot[categorical_var].astype(str)

    # Filter data
    count_df, unique_df = filter_and_count(
        data=data_to_plot,
        categorical_var=categorical_var,
        count_var=count_var,
        sig_var=sig_var,
        sig_threshold=sig_threshold
    )

    # Choose whether to plot all gene sets or only unique gene sets
    if unique:
        data_to_plot = unique_df
    else:
        data_to_plot = count_df

    if debug:
        print(f"Counts of motifs enriched: {data_to_plot}")

    # Example plot using Plotly (replace with your actual plotting code)
    fig = go.Figure(layout=dict(template='plotly'))
    fig = barplot(
        data=data_to_plot,
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
def update_num_enriched_traits_plot(
    selected_run, 
    debug=True,
    unique=False,
    categorical_var = "program_name",
    count_var = "trait_reported",
    sig_var = "adj_pval",
    sig_threshold = 0.05
):
    
    if debug:
        print(f"Selected run for trait enrichment plot: {selected_run}")
        print(f"Filtering trait enrichments with {sig_var} < {sig_threshold}")
        print(f"Counting unique {count_var} enriched for each {categorical_var}")

    # Assuming we want to plot something from the selected run
    data_to_plot = results['trait_enrichments'][selected_run]  # Adjust to your actual plot logic

    # Make sure x-axis is string
    data_to_plot[categorical_var] = data_to_plot[categorical_var].astype(str)

    # Filter data
    count_df, unique_df = filter_and_count(
        data=data_to_plot,
        categorical_var=categorical_var,
        count_var=count_var,
        sig_var=sig_var,
        sig_threshold=sig_threshold
    )

    # Choose whether to plot all gene sets or only unique gene sets
    if unique:
        data_to_plot = unique_df
    else:
        data_to_plot = count_df

    if debug:
        print(f"Counts of traits enriched: {data_to_plot}")

    # Example plot using Plotly (replace with your actual plotting code)
    fig = go.Figure(layout=dict(template='plotly'))
    fig = barplot(
        data=data_to_plot,
        x_column=categorical_var,
        y_column=count_var,
        title='',
        x_axis_title='Program',
        y_axis_title='Count',
        show_xaxis_labels=False
    )

    return fig


@callback(
    Output('num-perturbation-associations', 'figure'),
    [Input('run-selector', 'value')]
)
def update_num_perturbation_associations_plot(
    selected_run, 
    debug=True,
    unique=False,
    categorical_var = "program_name",
    count_var = "target_name",
    sig_var = "pval",
    sig_threshold = 0.25
):
        
    if debug:
        print(f"Selected run for perturbation association plot: {selected_run}")
        print(f"Filtering perturbation associations with {sig_var} < {sig_threshold}")
        print(f"Counting unique {count_var} enriched for each {categorical_var}")

    # Assuming we want to plot something from the selected run
    data_to_plot = results['perturbation_associations'][selected_run]

    # Make sure x-axis is string
    data_to_plot[categorical_var] = data_to_plot[categorical_var].astype(str)

    # Filter data
    count_df, unique_df = filter_and_count(
        data=data_to_plot,
        categorical_var=categorical_var,
        count_var=count_var,
        sig_var=sig_var,
        sig_threshold=sig_threshold
    )

    # Choose whether to plot all gene sets or only unique gene sets
    if unique:
        data_to_plot = unique_df
    else:
        data_to_plot = count_df

    if debug:
        print(f"Counts of perturbation associations: {data_to_plot}")

    # Example plot using Plotly (replace with your actual plotting code)
    fig = go.Figure(layout=dict(template='plotly'))
    fig = barplot(
        data=data_to_plot,
        x_column=categorical_var,
        y_column=count_var,
        title='',
        x_axis_title='Program',
        y_axis_title='Count',
        show_xaxis_labels=False
    )
    
    return fig


@callback(
    Output('phewas-binary-plot', 'figure'),
    Input('run-selector', 'value')
)
def update_phewas_binary_plot(selected_run):

    # Assuming we want to plot something from the selected run
    data_to_plot = results['trait_enrichments'][selected_run]

    fig = px.scatter(
        data_to_plot.query("trait_category != 'measurement'"),
        x='trait_reported',
        y='-log10(adj_pval)',
        color='trait_category',
        title="",
        hover_data=[
            "program_name", 
            "trait_reported", 
            "trait_category", 
            "adj_pval", 
            "genes", 
            "study_id", 
            "pmid"
        ]
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


@callback(
    Output('phewas-continuous-plot', 'figure'),
    Input('run-selector', 'value')
)
def update_phewas_continuous_plot(selected_run):

    # Assuming we want to plot something from the selected run
    data_to_plot = results['trait_enrichments'][selected_run]

    fig = px.scatter(
        data_to_plot.query("trait_category == 'measurement'"),
        x='trait_reported',
        y='-log10(adj_pval)',
        color='trait_category',
        title="",
        hover_data=[
            "program_name", 
            "trait_reported", 
            "trait_category", 
            "adj_pval", 
            "genes", 
            "study_id", 
            "pmid"
        ]
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
