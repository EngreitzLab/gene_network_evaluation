from dash import Input, Output
import plotly.express as px
from data_processing import count_unique
import pandas as pd

def register_callbacks(app, enrichment_gsea, enrichment_motif, enrichment_trait):
    @app.callback(
        Output('gsea-terms-table', 'data'),
        Input('programid-dropdown', 'value'),
        Input('qvalue-dropdown', 'value')
    )
    def update_gsea_table(selected_programids, selected_qvalue):
        filtered_df = enrichment_gsea[enrichment_gsea['qvalue'] <= selected_qvalue]
        if selected_programids:
            filtered_df = filtered_df[filtered_df['ProgramID'].isin(selected_programids)]
        return filtered_df.to_dict('records')

    @app.callback(
        Output('motif-terms-table', 'data'),
        Output('motif-stacked-barplot', 'figure'),
        Input('motif-programid-dropdown', 'value'),
        Input('ep-type-dropdown', 'value'),
        Input('motif-fdr-dropdown', 'value')
    )
    def update_motif_table(selected_programids, selected_ep_types, selected_fdr):
        filtered_df = enrichment_motif[enrichment_motif['FDR'] <= selected_fdr]
        if selected_programids:
            filtered_df = filtered_df[filtered_df['ProgramID'].isin(selected_programids)]
        if selected_ep_types:
            filtered_df = filtered_df[filtered_df['EPType'].isin(selected_ep_types)]
        motif_promoter_df = count_unique('ProgramID', 'TFMotif', filtered_df[filtered_df['EPType'] == 'Promoter'])
        motif_enhancer_df = count_unique('ProgramID', 'TFMotif', filtered_df[filtered_df['EPType'] == 'Enhancer'])
        motif_combined_df = pd.concat([
            motif_promoter_df.assign(Type='Promoter'),
            motif_enhancer_df.assign(Type='Enhancer')
        ])

        # Calculate total counts across promoters and enhancers
        program_totals = motif_combined_df.groupby('ProgramID')['TFMotif'].sum().reset_index(name='TotalMotifs')

        # Sort the programs by the total counts
        program_totals = program_totals.sort_values('TotalMotifs', ascending=False)

        # Merge the sorted program totals back to the combined dataframe
        motif_combined_df = motif_combined_df.merge(program_totals, on='ProgramID')

        fig = px.bar(
            motif_combined_df,
            x='ProgramID',
            y='TFMotif',
            color='Type',
            category_orders={'ProgramID': program_totals['ProgramID'].tolist()},
            barmode='stack',
            title='Motif Enrichment (Promoters and Enhancers)',
            template='simple_white'
        ).update_layout(xaxis_title='ProgramID', yaxis_title='Count')

        return filtered_df.to_dict('records'), fig

    @app.callback(
        Output('gwas-terms-table', 'data'),
        Output('gwas-barplot', 'figure'),
        Input('gwas-programid-dropdown', 'value'),
        Input('gwas-pvalue-dropdown', 'value')
    )
    def update_gwas_table(selected_programids, selected_pvalue):
        filtered_df = enrichment_trait[enrichment_trait['Adjusted P-value'] <= selected_pvalue]
        if selected_programids:
            filtered_df = filtered_df[filtered_df['ProgramID'].isin(selected_programids)]
        gwas_unique_df = count_unique('ProgramID', 'Term', filtered_df)
        fig = px.bar(
            gwas_unique_df.sort_values('Term', ascending=False),
            x='ProgramID',
            y='Term',
            title='GWAS Trait Enrichment',
            template='simple_white'
        ).update_layout(xaxis_title='ProgramID', yaxis_title='Count')
        return filtered_df.to_dict('records'), fig
