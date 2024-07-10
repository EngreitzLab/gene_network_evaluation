import pandas as pd
from dash import html, dcc
from dash.dash_table import DataTable
import plotly.express as px
import plotly.io as pio

# Set the default template to plotly_white
pio.templates.default = "plotly_white"

def create_layout(explained_variance, gsea_unique_df, enrichment_gsea, enrichment_motif, enrichment_trait):
    return html.Div([
        html.H1(children='Gene Program Evaluation Dashboard', style={'textAlign': 'center'}),
        html.Div(className="row", children=[
            html.Div(className="six columns", children=[
                dcc.Graph(
                    figure=px.scatter(
                        x=explained_variance.sort_values('VarianceExplained')['ProgramID'],
                        y=explained_variance.sort_values('VarianceExplained')['VarianceExplained'],
                        template='simple_white'
                    ).update_layout(title='Explained Variance by ProgramID', xaxis_title='ProgramID', yaxis_title='Variance Explained')
                ),
                dcc.Graph(
                    figure=px.bar(
                        gsea_unique_df.sort_values('ID', ascending=False),
                        x='ProgramID',
                        y='ID',
                        template='simple_white'
                    ).update_layout(title='Unique GSEA Terms per Program', xaxis_title='ProgramID', yaxis_title='Count')
                ),
            ]),
            html.Div(className="six columns", children=[
                html.H2("GSEA Terms Table"),
                html.Div([
                    html.Label('Filter by ProgramID'),
                    dcc.Dropdown(
                        id='programid-dropdown',
                        options=[{'label': i, 'value': i} for i in enrichment_gsea['ProgramID'].unique()],
                        multi=True
                    ),
                    html.Label('Filter by qvalue'),
                    dcc.Dropdown(
                        id='qvalue-dropdown',
                        options=[
                            {'label': '≤ 0.05', 'value': 0.05},
                            {'label': '≤ 0.01', 'value': 0.01},
                            {'label': '≤ 0.001', 'value': 0.001}
                        ],
                        multi=False,
                        value=0.05
                    )
                ]),
                DataTable(
                    id='gsea-terms-table',
                    columns=[{"name": i, "id": i} for i in enrichment_gsea.columns],
                    data=enrichment_gsea.to_dict('records'),
                    filter_action="native",
                    sort_action="native",
                    page_size=10
                )
            ])
        ]),
        html.Div(className="row", children=[
            html.Div(className="six columns", children=[
                html.H2("Motif Enrichment Table"),
                html.Div([
                    html.Label('Filter by ProgramID'),
                    dcc.Dropdown(
                        id='motif-programid-dropdown',
                        options=[{'label': i, 'value': i} for i in enrichment_motif['ProgramID'].unique()],
                        multi=True
                    ),
                    html.Label('Filter by EPType'),
                    dcc.Dropdown(
                        id='ep-type-dropdown',
                        options=[
                            {'label': 'Promoter', 'value': 'Promoter'},
                            {'label': 'Enhancer', 'value': 'Enhancer'}
                        ],
                        multi=True,
                        value=['Promoter', 'Enhancer']
                    ),
                    html.Label('Filter by FDR'),
                    dcc.Dropdown(
                        id='motif-fdr-dropdown',
                        options=[
                            {'label': '≤ 0.05', 'value': 0.05},
                            {'label': '≤ 0.01', 'value': 0.01},
                            {'label': '≤ 0.001', 'value': 0.001}
                        ],
                        multi=False,
                        value=0.05
                    )
                ]),
                dcc.Graph(id='motif-stacked-barplot'),
                DataTable(
                    id='motif-terms-table',
                    columns=[{"name": i, "id": i} for i in enrichment_motif.columns],
                    data=enrichment_motif.to_dict('records'),
                    filter_action="native",
                    sort_action="native",
                    page_size=10
                )
            ]),
            html.Div(className="six columns", children=[
                html.H2("GWAS Trait Enrichment Table"),
                html.Div([
                    html.Label('Filter by ProgramID'),
                    dcc.Dropdown(
                        id='gwas-programid-dropdown',
                        options=[{'label': i, 'value': i} for i in enrichment_trait['ProgramID'].unique()],
                        multi=True
                    ),
                    html.Label('Filter by Adjusted P-value'),
                    dcc.Dropdown(
                        id='gwas-pvalue-dropdown',
                        options=[
                            {'label': '≤ 0.05', 'value': 0.05},
                            {'label': '≤ 0.01', 'value': 0.01},
                            {'label': '≤ 0.001', 'value': 0.001}
                        ],
                        multi=False,
                        value=0.05
                    )
                ]),
                dcc.Graph(id='gwas-barplot'),
                DataTable(
                    id='gwas-terms-table',
                    columns=[{"name": i, "id": i} for i in enrichment_trait.columns],
                    data=enrichment_trait.to_dict('records'),
                    filter_action="native",
                    sort_action="native",
                    page_size=10
                )
            ])
        ])
    ])
