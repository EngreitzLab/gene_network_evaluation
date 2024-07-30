import dash
from dash import Dash, html, dcc, Input, Output
from data_processing import parse_mdata_summary, parse_mdata_cross_run, extract_total_unique_counts
from layout import create_scatter_layout, create_filtered_barplot_layout, create_filtered_stacked_barplot_layout
import dash_bootstrap_components as dbc
import mudata


# Load mudata
mdata = mudata.read_h5mu("../example_data/cNMF_evaluation_dashapp_data_small.h5mu")
summary_dict = parse_mdata_summary(mdata, verbose=False)
cross_run_dict = parse_mdata_cross_run(mdata, verbose=False)


# Create Dash app
app = Dash(__name__, use_pages=True, pages_folder="", external_stylesheets=[dbc.themes.SPACELAB])
app.title = "Gene Program Evaluation Dashboard v0.0.1"

# Summary page function: TODO move to pages/summary.py
def summary_page():
    return html.Div(id='landing_page', className='section', style={'width': '50%', 'display': 'inline-block', 'vertical-align': 'top'}, children=[
        html.H2("Input Data Summary"),
        html.P(f"{len(mdata.mod)} AnnData’s analyzed: {', '.join(mdata.mod.keys())}"),
        html.P(f"Number of cells in each AnnData: {summary_dict['n_cells']}"),
        html.P("Number of programs analyzed per modality:"),
        html.Ul([html.Li(f"{mod}: {summary_dict[mod]['n_programs']}") for mod in mdata.mod.keys()]),
    ])


# Cross run analysis: TODO move to pages/cross_run_analysis.py
totals = extract_total_unique_counts(cross_run_dict, summary_dict)
def cross_run_analysis():
    return html.Div(id='single_run_analysis', className='section', style={'width': '50%', 'display': 'inline-block', 'vertical-align': 'top'}, children=[
        html.H2("Cross Run Analysis"),
        html.Div([
            html.H2("Total Unique GSEA Terms enriched"),
            create_scatter_layout(
                data=totals,
                x_column='n_programs',
                y_column='enrichment_gsea',
                title='',
                x_axis_title='Number of Programs',
                y_axis_title='Unique GSEA Terms',
                id_suffix='gsea-cross_run',
            )
        ], style={'width': '100%', 'paddingBottom': '20px'}),
        html.Div([
            html.H2("Total Unique Motif Enrichment"),
            create_scatter_layout(
                data=totals,
                x_column='n_programs',
                y_column='enrichment_motif',
                title='',
                x_axis_title='Number of Programs',
                y_axis_title='Unique Motif Enrichment',
                id_suffix='motif-cross_run',
            )
        ], style={'width': '100%', 'paddingBottom': '20px'}),
        html.Div([
            html.H2("Total Unique GWAS Terms enriched"),
            create_scatter_layout(
                data=totals,
                x_column='n_programs',
                y_column='enrichment_gwas',
                title='',
                x_axis_title='Number of Programs',
                y_axis_title='Unique GWAS Terms',
                id_suffix='gwas-cross_run',
            )
        ], style={'width': '100%', 'paddingBottom': '20px'})
    ])

# Single run analysis: TODO move to pages/single_run_analysis.py
# TODO: make dropdown to select run --> currently failing to register filter barplot callback
# TODO: add in Narges's plots
curr_run = cross_run_dict["cNMF"]
def single_run_analysis():
    return html.Div(id='cross_run_analysis', className='section', style={'width': '50%', 'display': 'inline-block', 'vertical-align': 'top'}, children=[
        html.H2("Single Run Analysis"),
        html.P("Description: Analyze individual runs."),
        html.Div([
            html.H2("Explained Variance"),
            create_scatter_layout(
                data=curr_run["explained_variance"],
                x_column='ProgramID',
                y_column='VarianceExplained',
                title='',
                x_axis_title='Component',
                y_axis_title='Variance Explained (R²)',
                id_suffix='explained-variance'
            )
        ], style={'width': '100%', 'paddingBottom': '20px'}),

        html.Div([
            html.H2("Unique GSEA Terms enriched"),
            create_filtered_barplot_layout(
                app=app,
                data=curr_run["enrichment_gsea"],
                x_column='ProgramID',
                y_column='ID',
                filter_column='qvalue',
                starting_filter_value=0.05,
                title='',
                x_axis_title='Program',
                y_axis_title='Count',
                id_suffix='gsea-single_run',
                show_xaxis_labels=False
            )
        ], style={'width': '100%', 'paddingBottom': '20px'}),

        html.Div([
            html.H2("Motif Enrichment"),
            create_filtered_stacked_barplot_layout(
                app=app,
                data=curr_run["enrichment_motif"],
                x_column='ProgramID',
                y_column='TFMotif',
                stack_column='EPType',
                filter_column='FDR',
                starting_filter_value=0.05,
                title='',
                x_axis_title='Program',
                y_axis_title='Count',
                id_suffix='motif-single_run',
                show_xaxis_labels=False
            )
        ], style={'width': '100%', 'paddingBottom': '20px'}),

        html.Div([
            html.H2("Unique GWAS Terms enriched"),
            create_filtered_barplot_layout(
                app=app,
                data=curr_run["enrichment_gwas"],
                x_column='ProgramID',
                y_column='Term',
                filter_column='Adjusted_P-value',
                starting_filter_value=0.05,
                title='',
                x_axis_title='Program',
                y_axis_title='Count',
                id_suffix='gwas-single_run',
                show_xaxis_labels=False
            )
        ], style={'width': '100%', 'paddingBottom': '20px'})
    ])

# Program analysis: TODO move to pages/program_analysis.py
# TODO: add in single program plots and tables
def program_analysis():
    return html.Div(id='program_analysis', className='section', style={'width': '50%', 'display': 'inline-block', 'vertical-align': 'top'}, children=[
        html.H2("Program Analysis"),
        html.P("Description: Analyze individual programs."),
        html.P("This section will contain detailed analysis for a selected program.")
    ])


# Register pages
dash.register_page("Summary", path="/", layout=summary_page(), title="Summary",name="Summary")
dash.register_page("Cross Run Analysis", path="/cross-run", layout=cross_run_analysis(), title="Cross Run Analysis", name="Cross Run Analysis")
dash.register_page("Single Run Analysis", path="/single-run", layout=single_run_analysis(), title="Single Run Analysis", name="Single Run Analysis")
dash.register_page("Program Analysis", path="/program-analysis", layout=program_analysis(), title="Program Analysis", name="Program Analysis")
page_order = ["Summary", "Cross Run Analysis", "Single Run Analysis", "Program Analysis"]
dash.page_registry = {key: dash.page_registry[key] for key in page_order}

# Create sidebar for page navigation
sidebar = dbc.Nav(
    [
        dbc.NavLink(
            [
                html.Div(page["name"], className="sidebar-link"),
            ],
            href=page["relative_path"],
        )
        for page in dash.page_registry.values()
    ],
    vertical=True,
    pills=True,
    className="sidebar",
)
app.layout = app.layout = dbc.Container([
    dbc.Row([
        dbc.Col(html.Div("Gene Program Evaluation Dashboard v0.0.1",
                         style={'fontSize':50, 'textAlign':'center'}))
    ]),

    html.Hr(),

    dbc.Row(
        [
            dbc.Col(
                [
                    sidebar
                ], xs=4, sm=4, md=2, lg=2, xl=2, xxl=2),

            dbc.Col(
                [
                    dash.page_container
                ], xs=8, sm=8, md=10, lg=10, xl=10, xxl=10)
        ]
    )
], fluid=True)

# Run the app
if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0', port=8050)
