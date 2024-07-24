from dash import Input, Output, callback_context, dcc
from layout import create_barplot_layout, create_stacked_barplot_layout, filter_data, count_unique
import pandas as pd
import plotly.express as px



def register_callbacks(app, gsea_data, motif_data, gwas_data):
    @app.callback(
        Output('qvalue-barplot-container-gsea', 'children'),
        Output('FDR-stacked-barplot-container-motif', 'children'),
        Output('Adjusted P-value-barplot-container-gwas', 'children'),
        Input('qvalue-input-gsea', 'value'),
        Input('FDR-input-motif', 'value'),
        Input('Adjusted P-value-input-gwas', 'value')
    )
    def update_plots(gsea_qvalue, motif_fdr, gwas_pvalue):
        ctx = callback_context

        # Determine which input triggered the callback
        triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]

        if triggered_id == 'qvalue-input-gsea':
            filters = {'qvalue': gsea_qvalue}
            filtered_data = filter_data(gsea_data, filters)
            gsea_unique_df = count_unique(
                categorical_var='ProgramID',
                count_var='ID',
                dataframe=filtered_data
            )
            gsea_plot = create_barplot_layout(
                data=gsea_unique_df,
                x_column='ProgramID',
                y_column='ID',
                title='Unique GO Terms per Program',
                x_axis_title='ProgramID',
                y_axis_title='Count'
            )
        else:
            gsea_plot = dcc.Graph(figure={})

        if triggered_id == 'FDR-input-motif':
            filters = {'FDR': motif_fdr}
            filtered_data = filter_data(motif_data, filters)
            motif_promoter_df = count_unique('ProgramID', 'TFMotif', filtered_data[filtered_data['EPType'] == 'Promoter'])
            motif_enhancer_df = count_unique('ProgramID', 'TFMotif', filtered_data[filtered_data['EPType'] == 'Enhancer'])
            motif_combined_df = pd.concat([
                motif_promoter_df.assign(EPType='Promoter'),
                motif_enhancer_df.assign(EPType='Enhancer')
            ])
            motif_fig = px.bar(
                motif_combined_df.sort_values('TFMotif', ascending=False),
                x='ProgramID',
                y='TFMotif',
                color='EPType',
                barmode='stack',
                template='plotly_white'
            ).update_layout(
                title='Unique Enriched TF Motifs per Program',
                xaxis_title='ProgramID',
                yaxis_title='Count'
            )
            motif_plot = dcc.Graph(figure=motif_fig)
        else:
            motif_plot = dcc.Graph(figure={})

        if triggered_id == 'Adjusted P-value-input-gwas':
            filters = {'Adjusted P-value': gwas_pvalue}
            filtered_data = filter_data(gwas_data, filters)
            gwas_unique_df = count_unique(
                categorical_var='ProgramID',
                count_var='Term',
                dataframe=filtered_data
            )
            gwas_plot = create_barplot_layout(
                data=gwas_unique_df,
                x_column='ProgramID',
                y_column='Term',
                title='Unique GWAS Terms per Program',
                x_axis_title='ProgramID',
                y_axis_title='Count'
            )
        else:
            gwas_plot = dcc.Graph(figure={})

        return gsea_plot, motif_plot, gwas_plot
