import pandas as pd
import dash
from dash import Dash, html, dcc, Input, Output
from data_processing import load_data
from layout import create_scatter_layout, create_filtered_barplot_layout, create_filtered_stacked_barplot_layout
import dash_bootstrap_components as dbc


# Load data
explained_variance, enrichment_gsea, enrichment_motif, enrichment_gwas = load_data('../example_data/cNMF_evaluation_output.xlsx')
trait_phewas = pd.read_csv("../example_data/cNMF_enrichment_trait_processed.txt", sep='\t')

# Create Dash app
app = Dash(__name__, use_pages=True, pages_folder="", external_stylesheets=[dbc.themes.SPACELAB])
app.title = "Gene Program Evaluation Dashboard v0.0.1"


def summary_page():
    return html.Div(id='landing_page', className='section', style={'width': '50%', 'display': 'inline-block', 'vertical-align': 'top'}, children=[
        html.H2("Input Data Summary"),
        html.P("Number of AnnData’s analyzed: ..."),  # Placeholder
        html.P("Number of cells in each AnnData: ..."),  # Placeholder
        html.P("This section will be filled with actual data summary."),
    ])


def cross_run_analysis():
    return html.Div(id='cross_run_analysis', className='section', style={'width': '50%', 'display': 'inline-block', 'vertical-align': 'top'}, children=[
        html.H2("Cross Run Analysis"),
        html.P("Description: This section will compare different runs."),
        html.Div([
            html.H2("Explained Variance"),
            create_scatter_layout(
                data=explained_variance,
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
                data=enrichment_gsea,
                x_column='ProgramID',
                y_column='ID',
                filter_column='qvalue',
                starting_filter_value=0.05,
                title='',
                x_axis_title='Program',
                y_axis_title='Count',
                id_suffix='gsea-cross-run',
                show_xaxis_labels=False
            )
        ], style={'width': '100%', 'paddingBottom': '20px'}),

        html.Div([
            html.H2("Motif Enrichment"),
            create_filtered_stacked_barplot_layout(
                app=app,
                data=enrichment_motif,
                x_column='ProgramID',
                y_column='TFMotif',
                stack_column='EPType',
                filter_column='FDR',
                starting_filter_value=0.05,
                title='',
                x_axis_title='Program',
                y_axis_title='Count',
                id_suffix='motif-cross-run',
                show_xaxis_labels=False
            )
        ], style={'width': '100%', 'paddingBottom': '20px'}),

        html.Div([
            html.H2("Unique GWAS Terms enriched"),
            create_filtered_barplot_layout(
                app=app,
                data=enrichment_gwas,
                x_column='ProgramID',
                y_column='Term',
                filter_column='Adjusted P-value',
                starting_filter_value=0.05,
                title='',
                x_axis_title='Program',
                y_axis_title='Count',
                id_suffix='gwas-cross-run',
                show_xaxis_labels=False
            )
        ], style={'width': '100%', 'paddingBottom': '20px'})
    ])


def single_run_analysis():
    return html.Div(id='single_run_analysis', className='section', style={'width': '50%', 'display': 'inline-block', 'vertical-align': 'top'}, children=[
        html.H2("Single Run Analysis"),
        html.P("Description: Analyze individual runs."),
        html.P("This section will contain detailed analysis for a selected run.")
    ])


def program_analysis():
    return html.Div(id='program_analysis', className='section', style={'width': '50%', 'display': 'inline-block', 'vertical-align': 'top'}, children=[
        html.H2("Program Analysis"),
        html.P("Description: Analyze individual programs."),
        html.P("This section will contain detailed analysis for a selected program.")
    ])


dash.register_page(
    "Summary", 
    path="/", 
    layout=summary_page(), 
    name="Summary"
)
dash.register_page("Cross Run Analysis", path="/cross-run", layout=cross_run_analysis(), title="Cross Run Analysis")
dash.register_page("Single Run Analysis", path="/single-run", layout=single_run_analysis(), title="Single Run Analysis")
dash.register_page("Program Analysis", path="/program-analysis", layout=program_analysis(), title="Program Analysis")

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

if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0', port=8050)
