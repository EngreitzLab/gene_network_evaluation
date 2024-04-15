import pandas as pd
import numpy as np
import mudata
import plotly.express as px
from ipywidgets import interact, Dropdown, Output, VBox

def create_enrichment_plotting_df(mdata,
                       l2g_dataframe_path,
                       varm_key_name="FDR_GWAS",
                       uns_key_name='genesets_GWAS',
                       prog_key="cNMF"):
    
    # Load the GWAS dataframe
    l2g = pd.read_csv(l2g_dataframe_path, compression='gzip', low_memory=False)
    l2g_metadata = l2g[["trait_efos", "trait_reported", "trait_category"]].drop_duplicates()

    # Get P-values for each GWAS program
    prog_nam = mdata[prog_key].var.index 
    enrich_ps = pd.DataFrame(index=mdata[prog_key].uns[uns_key_name])

    for prog in prog_nam:
        # Extract p-values for the current program
        prog_p_values = mdata[prog_key][:, prog].varm[varm_key_name].flatten()
        
        # Add the p-values as a column in the DataFrame
        enrich_ps[prog] = pd.Series(prog_p_values, index=mdata[prog_key].uns[uns_key_name])

    #drop columns that are NA, meaning that enrichment could not be computed for that
    # Program x GWAS pair
    enrich_ps.dropna(axis=0, how='all', inplace=True)
    
    #merge with the L2G metadata to map trait_efos to trait_reported and trait_category
    enrich_ps = enrich_ps.reset_index(names='trait_efos').merge(l2g_metadata, on="trait_efos", how="left")
    
    #make the dataframe in long format
    enrich_ps = enrich_ps.melt(id_vars=["trait_efos", "trait_reported", "trait_category"], 
                               var_name="Program_Name", 
                               value_name="Enrichment_Pvals").drop_duplicates().sort_values(by=["trait_category", "Enrichment_Pvals"])
    
    #If the input P-value == 0, then replace it with the lowest non-zero P-value in the dataframe
    min_value = (enrich_ps.query("Enrichment_Pvals > 0")['Enrichment_Pvals']).min()
    
    #compute the -log(10) P-value and deal with edge-cases (e.g. P=0, P=1)
    enrich_ps.loc[enrich_ps["Enrichment_Pvals"] == 0, "Enrichment_Pvals"] = min_value #replace P=0 with min non-0 p-value
    enrich_ps['-log10(p-value)'] = abs(-1*np.log10(enrich_ps['Enrichment_Pvals']))
    
    enrich_ps = enrich_ps.query('trait_category != "measurement"')
    enrich_ps.reset_index(drop=True, inplace=True)
    
    return enrich_ps


def plot_interactive_phewas(data, x_column='trait_reported',
                            y_column='-log10(p-value)',
                            color_column='trait_category',
                            filter_column='Program_Name',
                            significance_threshold=0.05):
    
    # Get unique values for the filtering column
    filter_values = ['All'] + list(data[filter_column].unique())

    # Initialize output widget to display the plot
    output = Output()

    # Function to update plot based on dropdown selection
    def update_plot(selected_value):
        # Filter data based on selected value
        if selected_value == "All":
            filtered_data = data.copy()  # No selection, show all data
        else:
            filtered_data = data[data[filter_column] == selected_value]

        # Create the plot
        fig = px.scatter(filtered_data, x=x_column, y=y_column, color=color_column,
                         title='Cell Program x OpenTargets GWAS L2G Enrichment',
                         hover_data=[filter_column, 'trait_efos', x_column, color_column, 'Enrichment_Pvals'])

        # Customize layout
        fig.update_layout(
            xaxis_title=x_column,
            yaxis_title=y_column,
            yaxis=dict(tickformat=".1f"),
            width=1000,  # Adjust width as needed
            height=500,  # Adjust height as needed,
            xaxis_tickfont=dict(size=4)
        )

        # Add horizontal dashed line for significance threshold
        fig.add_hline(y=-np.log10(significance_threshold), line_dash="dash",
                      annotation_text=f'Significance Threshold ({significance_threshold})', annotation_position="top right")

        # Clear previous plot and display the new one
        with output:
            output.clear_output(wait=True)
            fig.show()

    # Create dropdown widget
    dropdown = Dropdown(options=filter_values, description=f"{filter_column}:")

    # Define function to handle dropdown value change
    def on_change(change):
        if change['type'] == 'change' and change['name'] == 'value':
            update_plot(change['new'])

    # Link dropdown change to function
    dropdown.observe(on_change)

    # Display dropdown widget and initial plot
    display(VBox([dropdown, output]))