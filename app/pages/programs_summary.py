import os
import dash
import pickle
from dash import html, dcc, callback, Input, Output
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import pandas as pd
from plot import fig_to_uri, scatterplot, barplot, plot_topic_trait_relationship_heatmap
from utils import filter_and_count
import diskcache
import matplotlib.pyplot as plt


# Register the page
dash.register_page(__name__, order=1)

# Get the app and cache
app = dash.get_app()
cache = diskcache.Cache("./.cache")
results = cache.get("results")

# Get the first data type and its corresponding second-level keys as default
default_data_type = list(results.keys())[0]
default_run = list(results[default_data_type].keys())[0]

# Get the columns that don't have the selected dim reduction prefix
categorical_keys = results["categorical_keys"]
covariate_keys = categorical_keys
default_covariate = categorical_keys[0] if categorical_keys else None

# Create layout
layout = dbc.Container([

    # Header
    html.H1("Programs Summary", className="mb-4"),
    html.P(
        "This page is designed to provide an overview of the evaluations run on a given set of programs. "
        "It is broken up into three main sections: Goodness of Fit, Covariate Association, and Trait Enrichment. "
        "Goodness of Fit gives an overview of the number of enrichments and associations detected across programs at a user-defined threshold, "
        "and can be useful for determining if the value of k selected is appropriate for the data. "
        "Looking at covariate associations of programs can help identify if the programs are associated with technical or biological sources of variation. "
        "Trait enrichment can help identify if the programs are associated with any traits of interest."
    ),
           
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
            html.P("This tab provides an overview the number of enrichments and associations detected across programs at a user define p-value threshold."),
            
            # Start with a significance threshold input
            dbc.Row([
                dbc.Col([
                    html.H3("Select threshold"),
                    html.P("Select the significance threshold to use for filtering enrichments and associations."),
                    dcc.Input(
                        id='terms-threshold',
                        type='number',
                        value=0.05,
                        min=0,
                        max=1,
                        placeholder="Enter significance threshold"
                    ),
                ], width=6),
            ], className="mb-4"),
            
            # Explained Variance Per Program and Number of Genesets Enriched Per Program
            dbc.Row([

                # Explained Variance Per Program
                dbc.Col([
                    html.H3("Explained Variance Per Program"),
                    dcc.Graph(id='explained-variance-plot'),
                ], width=6),

                # Number of Genesets Enriched Per Program
                dbc.Col([
                    html.H3("Number of Genesets Enriched Per Program"),
                    dcc.Graph(id='num-enriched-genesets'),
                ], width=6),
            ], className="mb-4"),

            # Number of Motifs Enriched Per Program and Number of Traits Enriched Per Program
            dbc.Row([

                # Number of Motifs Enriched Per Program
                dbc.Col([
                    html.H3("Number of Motifs Enriched Per Program"),
                    dcc.Graph(id='num-enriched-motifs'),
                ], width=6),

                # Number of Traits Enriched Per Program
                dbc.Col([
                    html.H3("Number of Traits Enriched Per Program"),
                    dcc.Graph(id='num-enriched-traits'),
                ], width=6),
            ], className="mb-4"),

            # Number of Perturbation Associations Per Program
            dbc.Row([
                # Number of Perturbation Associations Per Program
                dbc.Col([
                    html.H3("Number of Perturbation Associations Per Program", className="mt-4"),
                    dcc.Graph(id='num-perturbation-associations'),
                ], width=6),
            ], className="mb-4"),
        ]),

        # 2. Tab for Covariate Association
        dcc.Tab(label='Categorical Association', children=[
            html.H2("Categorical Association", className="mt-3 mb-3"),
            html.P("This tab provides an overview of the categorical associations between programs and a selected covariate."),
            
            # Covariate dropdown
            dbc.Row([
                dbc.Col([
                    html.Label("Select Covariate"),
                    dcc.Dropdown(
                        id='covariate-selector',
                        options=[{'label': covariate, 'value': covariate} for covariate in covariate_keys],
                        value=default_covariate
                    )
                ], width=6),
            ], className="mb-4"),

            # Program trait correlation heatmap
            dbc.Row([
                dbc.Col([
                    html.H3("Program trait correlation"),
                    html.P("This heatmap shows the correlations between programs and a passed in categorical trait."),
                    html.Div([html.Img(id='program-trait-heatmap', className="mt-3")])
                ], width=12),
            ], className="mb-4"),

            # Volcano plot for enrichment with respect to covariate
            dbc.Row([
                dbc.Col([
                    html.H3("Categorical association testing"),
                    html.P("A categorical variable walks into a bar..."),
                    dcc.Graph(id='categorical-assoc-volcano-plot'),
                ], width=12),
            ], className="mb-4"),
        ]),

        # 3. Tab for Trait Enrichment
        dcc.Tab(label='Trait Enrichment', children=[
            html.H2("Trait Enrichment", className="mt-3 mb-3"),
            html.P("This tab provides two PheWAS style plots for investigating the enrichment of gene programs with binary and continuous traits. "
                   "The x-axis represents the trait reported, the y-axis represents the -log10 adjusted p-value, and the color represents the trait category."),

            # Binary trait enrichment PheWAS plot
            dbc.Row([
                dbc.Col([
                    html.H3("Binary trait enrichment PheWAS plot"),
                    dcc.Graph(id='phewas-binary-plot'),
                ], width=12),
            ], className="mb-4"),

            # Continuous trait enrichment PheWAS plot
            dbc.Row([
                dbc.Col([
                    html.H3("Continuous trait enrichment PheWAS plot"),
                    dcc.Graph(id='phewas-continuous-plot'),
                ], width=12),
            ], className="mb-4"),
        ]),
    
    ])

], fluid=True, className="p-4")


# Callback for explained variance
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
    
    # make all the colors blue
    colors = ["blue"] * len(data_to_plot)

    if debug:
        print(f"Explained variance ratios: {data_to_plot}")

    # Example plot using Plotly (replace with your actual plotting code)
    fig = scatterplot(
        data=data_to_plot,
        x_column='program_name',
        y_column='explained_variance_ratio',
        title='',
        x_axis_title='Component',
        y_axis_title='Variance Explained (R²)',
        colors=colors,
        size=8
    )
    return fig


# Callback for number of enriched gene sets
@callback(
    Output('num-enriched-genesets', 'figure'),
    [Input('run-selector', 'value'), Input('terms-threshold', 'value')]
)
def update_num_enriched_genesets_plot(
    selected_run, 
    sig_threshold,
    debug=False,
    unique=False,
    categorical_var = "program_name",
    count_var = "term",
    sig_var = "adj_pval",
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


# Callback for number of enriched motifs
@callback(
    Output('num-enriched-motifs', 'figure'),
    [Input('run-selector', 'value'), Input('terms-threshold', 'value')]
)
def update_num_enriched_motifs_plot(
    selected_run, 
    sig_threshold,
    debug=False,
    unique=False,
    categorical_var = "program_name",
    count_var = "motif",
    sig_var = "pval",
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


# Callback for number of enriched traits
@callback(
    Output('num-enriched-traits', 'figure'),
    [Input('run-selector', 'value'), Input('terms-threshold', 'value')]
)
def update_num_enriched_traits_plot(
    selected_run,
    sig_threshold,
    debug=False,
    unique=False,
    categorical_var = "program_name",
    count_var = "trait_reported",
    sig_var = "adj_pval",
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


# Callback for number of perturbation associations
@callback(
    Output('num-perturbation-associations', 'figure'),
    [Input('run-selector', 'value'), Input('terms-threshold', 'value')]
)
def update_num_perturbation_associations_plot(
    selected_run, 
    sig_threshold,
    debug=False,
    unique=False,
    categorical_var = "program_name",
    count_var = "target_name",
    sig_var = "pval",
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


# Callback for program trait correlation heatmap
@callback(
    Output('program-trait-heatmap', 'src'),
    [Input('run-selector', 'value'), Input('covariate-selector', 'value')]
)
def update_program_trait_heatmap(
    selected_run, 
    selected_covariate,
    debug=True
):
        
    if debug:
        print(f"Selected run for program trait heatmap: {selected_run}")
        print(f"Selected covariate for program trait heatmap: {selected_covariate}")

    # Assuming we want to plot something from the selected run
    obs = results["obs"]
    unique_values = obs[selected_covariate].unique()

    # Make sure x-axis is string
    curr_obs_membership = results['obs_memberships'][selected_run]
    n_programs = curr_obs_membership.shape[1]

    # Create a figure
    fig, ax = plt.subplots(figsize=(n_programs * 0.5, len(unique_values) * 1))

    # Example plot using Plotly (replace with your actual plotting code)
    plot_topic_trait_relationship_heatmap(
        ax,  
        curr_obs_membership,
        obs,
        covariates=[selected_covariate],
    )

    # Save the plot to a temporary file
    fig = fig_to_uri(fig)

    return fig


# Callback for categorical association volcano plot
@callback(
    Output('categorical-assoc-volcano-plot', 'figure'),
    [Input('run-selector', 'value'), Input('covariate-selector', 'value')]
)
def update_categorical_association_volcano_plot(
    selected_run, 
    selected_covariate,
    debug=True,
    sig_var = "_kruskall_wallis_pval",
    sig_threshold = 0.05,
    effect_var = "_kruskall_wallis_stat",
    effect_threshold = 0.5
    
):

    if debug:
        print(f"Selected run for categorical association plot: {selected_run}")
        print(f"Selected covariate for categorical association plot: {selected_covariate}")    

    # Check if the selected covariate is in the categorical associations, if it isn't return an empty figure with title
    # "Categorical association evaluation not run on selected covariate, no data to display"
    if f"{selected_covariate}{sig_var}" not in results['categorical_associations'][selected_run].columns:
        print(f"Selected covariate {selected_covariate} not found in categorical associations")
        fig = go.Figure(layout=dict(template='plotly'))
        fig.update_layout(title="Categorical association evaluation not run on selected covariate, no data to display")
        return fig

    # Assuming we want to plot something from the selected run
    data_to_plot = results['categorical_associations'][selected_run]
    data_to_plot["program_name"] = results["obs_memberships"][selected_run].columns

    # Filter data
    effect_var = f"{selected_covariate}{effect_var}"
    sig_var = f"{selected_covariate}{sig_var}"

    # -log10 transform the p-values, add a small value to avoid log(0)
    data_to_plot[f"-log10({sig_var})"] = -np.log10(data_to_plot[sig_var] + 1e-10)

    if debug:
        print(f"Data to plot: {data_to_plot}")

    # Plot a volcano plot with the effect size on the x-axis and the -log10(p-value) on the y-axis
    # Highlight those points that are significant and have a large effect size
    fig = px.scatter(
        data_to_plot,
        x=effect_var,
        y=f"-log10({sig_var})",
        color=sig_var,
        title="",
        hover_data=[
            "program_name",
            effect_var,
            sig_var
        ]
    )

    # Customize layout
    fig.update_layout(
        xaxis_title='Effect Size',
        yaxis_title='-log10(p-value)',
        yaxis=dict(tickformat=".1f"),
        width=1200,
        height=600,
        xaxis_tickfont=dict(size=4),
        plot_bgcolor='rgba(0,0,0,0)',  # Transparent background
    )

    # Add horizontal dashed line for significance threshold
    fig.add_hline(
        y=-np.log10(sig_threshold),
        line_dash="dash",
        annotation_text='Significance Threshold (0.05)',
        annotation_position="top right"
    )

    # Add vertical dashed line for effect size threshold
    fig.add_vline(
        x=effect_threshold,
        line_dash="dash",
        annotation_text='Effect Size Threshold',
        annotation_position="top right"
    )

    return fig

    
# Callback for phewas binary plot
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
        width=1200,
        height=600,
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
