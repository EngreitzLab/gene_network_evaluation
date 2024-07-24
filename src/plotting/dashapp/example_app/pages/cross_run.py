import dash
from dash import html, dcc, callback, Input, Output
from layout import create_scatter_layout, create_filtered_barplot_layout, create_filtered_stacked_barplot_layout

def cross_run_analysis(unique_id):
    return html.Div(id=f'cross_run_analysis_{unique_id}', className='section', style={'width': '50%', 'display': 'inline-block', 'vertical-align': 'top'}, children=[
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
                id_suffix=f'explained-variance-{unique_id}'
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
                x_axis_title='ProgramID',
                y_axis_title='Count',
                id_suffix=f'gsea-cross-run-{unique_id}',
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
                x_axis_title='ProgramID',
                y_axis_title='Count',
                id_suffix=f'motif-cross-run-{unique_id}',
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
                x_axis_title='ProgramID',
                y_axis_title='Count',
                id_suffix=f'gwas-cross-run-{unique_id}',
                show_xaxis_labels=False
            )
        ], style={'width': '100%', 'paddingBottom': '20px'})
    ])