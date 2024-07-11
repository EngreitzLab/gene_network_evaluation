# In callbacks.py
from dash import Input, Output, State
import plotly.express as px
from data_processing import count_unique
import pandas as pd
import numpy as np

def register_callbacks(app, enrichment_gsea, enrichment_motif, enrichment_trait, phewas_data):
    @app.callback(
        Output('gsea-terms-table', 'data'),
        Output('gsea-barplot', 'figure'),
        Input('programid-dropdown', 'value'),
        Input('qvalue-input', 'value')
    )
    def update_gsea_table(selected_programids, selected_qvalue):
        filtered_df = enrichment_gsea[enrichment_gsea['qvalue'] <= selected_qvalue]
        if selected_programids:
            filtered_df = filtered_df[filtered_df['ProgramID'].isin(selected_programids)]
        gsea_unique_df = count_unique(
            categorical_var='ProgramID',
            count_var='ID',
            dataframe=filtered_df
        )
        fig = px.bar(
            gsea_unique_df.sort_values('ID', ascending=False),
            x='ProgramID',
            y='ID',
            template="plotly_white"
        ).update_layout(title='Unique GSEA Terms per Program', xaxis_title='ProgramID', yaxis_title='Count')
        return filtered_df.to_dict('records'), fig

    @app.callback(
        Output('motif-terms-table', 'data'),
        Output('motif-stacked-barplot', 'figure'),
        Input('motif-programid-dropdown', 'value'),
        Input('ep-type-dropdown', 'value'),
        Input('motif-fdr-input', 'value')
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

        # Merge the total counts back to the combined dataframe
        motif_combined_df = motif_combined_df.merge(program_totals, on='ProgramID')

        # Sort the programs by the total counts
        program_totals = program_totals.sort_values('TotalMotifs', ascending=False)
        motif_combined_df = motif_combined_df.sort_values('TotalMotifs', ascending=False)

        fig = px.bar(
            motif_combined_df,
            x='ProgramID',
            y='TFMotif',
            color='Type',
            category_orders={'ProgramID': program_totals['ProgramID'].tolist()},
            barmode='stack',
            template="plotly_white"
        ).update_layout(title='Motif Enrichment (Promoters and Enhancers)', xaxis_title='ProgramID', yaxis_title='Count')
        return filtered_df.to_dict('records'), fig

    @app.callback(
        Output('gwas-terms-table', 'data'),
        Output('gwas-barplot', 'figure'),
        Input('gwas-programid-dropdown', 'value'),
        Input('gwas-pvalue-input', 'value')
    )
    def update_gwas_table(selected_programids, selected_pvalue):
        filtered_df = enrichment_trait[enrichment_trait['Adjusted P-value'] <= selected_pvalue]
        if selected_programids:
            filtered_df = filtered_df[filtered_df['ProgramID'].isin(selected_programids)]
        gwas_unique_df = count_unique(
            categorical_var='ProgramID',
            count_var='Term',
            dataframe=filtered_df
        )
        fig = px.bar(
            gwas_unique_df.sort_values('Term', ascending=False),
            x='ProgramID',
            y='Term',
            template="plotly_white"
        ).update_layout(title='Unique GWAS Traits per Program', xaxis_title='ProgramID', yaxis_title='Count')
        return filtered_df.to_dict('records'), fig

    @app.callback(
        Output('phewas-plot-binary', 'figure'),
        Input('phewas-program-dropdown-binary', 'value'),
        State('phewas-data', 'data')
    )
    def update_phewas_plot_binary(selected_program, data):
        data = pd.DataFrame(data)
        filtered_data = data.query("trait_category != 'measurement'")
        if selected_program != "All":
            filtered_data = filtered_data[filtered_data['program_name'] == selected_program]

        fig = px.scatter(
            filtered_data,
            x='trait_reported',
            y='-log10(p-value)',
            color='trait_category',
            title="Endothelial Cell Programs x GWAS Binary Outcome Enrichments",
            hover_data=["program_name", "trait_reported", "trait_category", "P-value", "Genes", "study_id", "pmid"]
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

    @app.callback(
        Output('phewas-plot-continuous', 'figure'),
        Input('phewas-program-dropdown-continuous', 'value'),
        State('phewas-data', 'data')
    )
    def update_phewas_plot_continuous(selected_program, data):
        data = pd.DataFrame(data)
        filtered_data = data.query("trait_category == 'measurement'")
        if selected_program != "All":
            filtered_data = filtered_data[filtered_data['program_name'] == selected_program]

        fig = px.scatter(
            filtered_data,
            x='trait_reported',
            y='-log10(p-value)',
            color='trait_category',
            title="Endothelial Cell Programs x GWAS Continuous Outcome Enrichments",
            hover_data=["program_name", "trait_reported", "trait_category", "P-value", "Genes", "study_id", "pmid"]
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

