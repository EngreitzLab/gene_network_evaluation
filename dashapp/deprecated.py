
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





def create_explained_variance_layout(
    explained_variance: pd.DataFrame,
    cumulative: bool = False
):
    """Create the layout for the explained variance plot.
    
    Parameters
    ----------
    explained_variance : pd.DataFrame
        DataFrame containing the explained variance data. Expected columns are 'ProgramID' and 'VarianceExplained'.

    Returns
    -------
    html.Div
        A Div containing the explained variance plot.
    """
    if cumulative:
        explained_variance['VarianceExplained'] = explained_variance['VarianceExplained'].cumsum()
    return html.Div([
        dcc.Graph(
            figure=px.scatter(
                x=explained_variance.sort_values('VarianceExplained', ascending=cumulative)['ProgramID'],
                y=explained_variance.sort_values('VarianceExplained', ascending=cumulative)['VarianceExplained'],
                template='plotly_white'
            ).update_layout(title='Explained Variance by ProgramID', xaxis_title='Component', yaxis_title='Variance Explained (R²)')
        )
    ])


def create_gsea_layout(
    enrichment_gsea: pd.DataFrame,
    starting_fdr: float = 0.05,
    cumulative: bool = False
):
    """Create the layout for the GSEA table and plot.

    Parameters
    ----------
    enrichment_gsea : pd.DataFrame
        DataFrame containing the GSEA enrichment data. Expected columns are 'ProgramID', 'ID', 'qvalue'.

    Returns
    -------
    html.Div
        A Div containing the GSEA table and plot.
    """
    # Preprocess data
    gsea_unique_df = count_unique(
        categorical_var='ProgramID',
        count_var='ID',
        dataframe=enrichment_gsea.loc[enrichment_gsea['qvalue'] <= starting_fdr]
    )
    if cumulative:
        gsea_unique_df['ID'] = gsea_unique_df['ID'].fillna(gsea_unique_df['ID'].max())
    else:
        gsea_unique_df['ID'] = gsea_unique_df['ID'].fillna(0)

    return html.Div([
        html.H2("Enriched GSEA Terms"),
        html.Div([
            html.Label('Filter by program'),
            dcc.Dropdown(
                id='programid-dropdown',
                options=[{'label': i, 'value': i} for i in enrichment_gsea['ProgramID'].unique()],
                multi=True
            ),
            html.Label('Filter by qvalue'),
            # Allow filtering by ANY user input q-value
            dcc.Input(
                id='qvalue-input',
                type='number',
                value=starting_fdr,
                step=0.001
            )
        ]),
        DataTable(
            id='gsea-terms-table',
            columns=[{"name": i, "id": i} for i in enrichment_gsea.columns],
            data=enrichment_gsea.to_dict('records'),
            filter_action="native",
            sort_action="native",
            page_size=10
        ),
        dcc.Graph(
            id='gsea-barplot',
            figure=px.bar(
                gsea_unique_df.sort_values('ID', ascending=False),
                x='ProgramID',
                y='ID',
                template="plotly_white"
            ).update_layout(title='Unique GSEA Terms per Program', xaxis_title='ProgramID', yaxis_title='Count')
        )
    ])


def create_motif_enrichment_layout(
    enrichment_motif: pd.DataFrame,
    starting_fdr: float = 0.05,
    cumulative: bool = False
):
    """Create the layout for the motif enrichment table and plot.

    Parameters
    ----------
    enrichment_motif : pd.DataFrame
        DataFrame containing the motif enrichment data. Expected columns are 'ProgramID', 'TFMotif', 'FDR', 'EPType'.

    Returns
    -------
    html.Div
        A Div containing the motif enrichment table and plot.
    """
    # Preprocess data
    filtered_df = enrichment_motif[enrichment_motif['FDR'] <= starting_fdr]
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

    return html.Div([
        html.H2("Motif Enrichment"),
        html.Div([
            html.Label('Filter by program'),
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
            dcc.Input(
                id='motif-fdr-input',
                type='number',
                value=starting_fdr,
                step=0.001
            )
        ]),
        DataTable(
            id='motif-terms-table',
            columns=[{"name": i, "id": i} for i in enrichment_motif.columns],
            data=enrichment_motif.to_dict('records'),
            filter_action="native",
            sort_action="native",
            page_size=10
        ),
        dcc.Graph(
            id='motif-stacked-barplot',
            figure=px.bar(
                motif_combined_df,
                x='ProgramID',
                y='TFMotif',
                color='Type',
                category_orders={'ProgramID': program_totals['ProgramID'].tolist()},
                barmode='stack',
                template="plotly_white"
            ).update_layout(title='Motif Enrichment (Promoters and Enhancers)', xaxis_title='ProgramID', yaxis_title='Count')
        )
    ])


def create_gwas_layout(
    enrichment_trait: pd.DataFrame,
    starting_pvalue: float = 0.05,
    cumulative: bool = False
):
    """Create the layout for the GWAS trait enrichment table and plot.

    Parameters
    ----------
    enrichment_trait : pd.DataFrame
        DataFrame containing the GWAS trait enrichment data. Expected columns are 'ProgramID', 'Term', 'Adjusted P-value'.

    Returns
    -------
    html.Div
        A Div containing the GWAS trait enrichment table and plot.
    """
    # Preprocess data
    gwas_unique_df = count_unique(
        categorical_var='ProgramID',
        count_var='Term',
        dataframe=enrichment_trait.loc[enrichment_trait['Adjusted P-value'] <= starting_pvalue]
    )
    if cumulative:
        gwas_unique_df['Term'] = gwas_unique_df['Term'].fillna(gwas_unique_df['Term'].max())
    else:
        gwas_unique_df['Term'] = gwas_unique_df['Term'].fillna(0)

    return html.Div([
        html.H2("Enriched GWAS Traits"),
        html.Div([
            html.Label('Filter by program'),
            dcc.Dropdown(
                id='gwas-programid-dropdown',
                options=[{'label': i, 'value': i} for i in enrichment_trait['ProgramID'].unique()],
                multi=True
            ),
            html.Label('Filter by adjusted p-value'),
            # Allow filtering by ANY user input p-value
            dcc.Input(
                id='gwas-pvalue-input',
                type='number',
                value=starting_pvalue,
                step=0.001
            )
        ]),
        DataTable(
            id='gwas-terms-table',
            columns=[{"name": i, "id": i} for i in enrichment_trait.columns],
            data=enrichment_trait.to_dict('records'),
            filter_action="native",
            sort_action="native",
            page_size=10
        ),
        dcc.Graph(
            id='gwas-barplot',
            figure=px.bar(
                gwas_unique_df.sort_values('Term', ascending=False),
                x='ProgramID',
                y='Term',
                template="plotly_white"
            ).update_layout(title='Unique GWAS Traits per Program', xaxis_title='ProgramID', yaxis_title='Count')
        )
    ])


def create_phewas_layout(
    phewas_data: pd.DataFrame,
    title: str = "PheWAS Plot",
    query_string: str = None,
    id_suffix: str = ""
):
    """Create the layout for the PheWAS plot.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame containing the PheWAS data.
    title : str
        Title of the plot.
    query_string : str, optional
        Query string to filter the data.
    id_suffix : str, optional
        Suffix to add to the id of the dropdown and plot.

    Returns
    -------
    html.Div
        A Div containing the PheWAS plot and dropdown.
    """
    if query_string:
        phewas_data = phewas_data.query(query_string)

    filter_column = 'program_name'
    filter_values = ['All'] + sorted(phewas_data[filter_column].unique())

    return html.Div([
        html.H2(title),
        html.Label('Filter by program'),
        dcc.Dropdown(
            id=f'phewas-program-dropdown-{id_suffix}',
            options=[{'label': i, 'value': i} for i in filter_values],
            value='All',
            multi=False
        ),
        dcc.Graph(id=f'phewas-plot-{id_suffix}')
    ])


def create_layout(
    explained_variance: pd.DataFrame,
    enrichment_gsea: pd.DataFrame,
    enrichment_motif: pd.DataFrame,
    enrichment_trait: pd.DataFrame,
    phewas_data: pd.DataFrame
):
    """Create the layout for the Dash app.

    Parameters
    ----------
    explained_variance : pd.DataFrame
        DataFrame containing the explained variance data. Expected columns are 'ProgramID' and 'VarianceExplained'.
    enrichment_gsea : pd.DataFrame
        DataFrame containing the GSEA enrichment data. Expected columns are 'ProgramID', 'ID', 'qvalue'.
    enrichment_motif : pd.DataFrame
        DataFrame containing the motif enrichment data. Expected columns are 'ProgramID', 'TFMotif', 'FDR', 'EPType'.
    enrichment_trait : pd.DataFrame
        DataFrame containing the GWAS trait enrichment data. Expected columns are 'ProgramID', 'Term', 'Adjusted P-value'.
    phewas_data : pd.DataFrame
        DataFrame containing the PheWAS data. Expected columns are 'program_name', 'trait_reported', 'trait_category', 'P-value', 'Genes', 'study_id

    Returns
    -------
    html.Div
        A Div containing the layout for the Dash app.
    """
    return html.Div([
        dcc.Store(id='phewas-data', data=phewas_data.to_dict('records')),
        html.H1("Gene Program Evaluation Dashboard", style={'textAlign': 'center'}),
        
        # Explained Variance Section
        html.Div(id='explained-variance-section', className='section', children=[
            html.H2("Explained Variance"),
            create_explained_variance_layout(explained_variance)
        ]),
        
        # Gene Set Enrichment Section
        html.Div(id='gsea-section', className='section', children=[
            html.H2("Gene Set Enrichment"),
            create_gsea_layout(enrichment_gsea)
        ]),
        
        # Motif Enrichment Section
        html.Div(id='motif-enrichment-section', className='section', children=[
            html.H2("Motif Enrichment"),
            create_motif_enrichment_layout(enrichment_motif)
        ]),
        
        # GWAS Enrichment Section
        html.Div(id='gwas-enrichment-section', className='section', children=[
            html.H2("GWAS Enrichment"),
            create_gwas_layout(enrichment_trait),
            create_phewas_layout(phewas_data, title="Endothelial Cell Programs x GWAS Binary Outcome Enrichments", query_string="trait_category != 'measurement'", id_suffix="binary"),
            create_phewas_layout(phewas_data, title="Endothelial Cell Programs x GWAS Continuous Outcome Enrichments", query_string="trait_category == 'measurement'", id_suffix="continuous")
        ]),
        
        # Placeholder Sections for Future Content
        html.Div(id='future-content-section', className='section', children=[
            html.H2("Future Content"),
            
            html.Div(id='preprocessing-section', className='section', children=[
                html.H3("Preprocessing"),
                html.P("Placeholder for preprocessing content.", className='placeholder')
            ]),
            
            html.Div(id='covariate-associations-section', className='section', children=[
                html.H3("Covariate Associations"),
                html.P("Placeholder for covariate associations content.", className='placeholder')
            ]),
            
            html.Div(id='program-loadings-section', className='section', children=[
                html.H3("Program Loadings"),
                html.P("Placeholder for program loadings content.", className='placeholder')
            ]),
            
            html.Div(id='regulators-section', className='section', children=[
                html.H3("Regulators"),
                html.P("Placeholder for regulators content.", className='placeholder')
            ]),
            
            html.Div(id='program-curation-section', className='section', children=[
                html.H3("Program Curation"),
                html.P("Placeholder for program curation content.", className='placeholder')
            ]),
        ])
    ])
