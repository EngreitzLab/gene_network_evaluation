import pandas as pd
import numpy as np
import mudata
import plotly.express as px
from ipywidgets import interact, Dropdown, Output, VBox

import pandas as pd
import numpy as np

def process_enrichment_data(enrich_res,
                            metadata,
                            pval_col="P-value",
                            enrich_geneset_id_col="Term",
                            metadata_geneset_id_col="trait_efos",
                            color_category_col="trait_category",
                            program_name_col="program_name",
                            annotation_cols=["trait_reported", "Genes", "study_id", "pmid"]):

    # Read in enrichment results
    if isinstance(enrich_res, str):
        enrich_df = pd.read_csv(enrich_res)
    elif isinstance(enrich_res, pd.DataFrame):
        enrich_df = enrich_res
    else:
        raise ValueError("enrich_res must be either a pandas DataFrame or a file path to a CSV file.")

    if isinstance(metadata, str):
        metadata_df = pd.read_csv(metadata, compression='gzip', low_memory=False)
    elif isinstance(metadata, pd.DataFrame):
        metadata_df = metadata
    else:
        raise ValueError("metadata must be either a pandas DataFrame or a file path to a CSV file.")

    # Join the enrichment results and the metadata
    enrich_ps = enrich_df.merge(metadata_df, left_on=enrich_geneset_id_col, right_on=metadata_geneset_id_col, how="left")

    # Only keep the relevant columns
    keep_cols = list([enrich_geneset_id_col, pval_col, metadata_geneset_id_col, color_category_col, program_name_col] + annotation_cols)
    enrich_ps = enrich_ps[keep_cols]

    # Sort by P-value
    enrich_ps = enrich_ps.drop_duplicates().sort_values(by=[color_category_col, pval_col])

    # If the input P-value == 0, then replace it with the lowest non-zero P-value in the dataframe
    min_value = enrich_ps.query(f"`{pval_col}` > 0")[pval_col].min()

    # Compute the -log(10) P-value and deal with edge-cases (e.g. P=0, P=1)
    enrich_ps.loc[enrich_ps[pval_col] == 0, pval_col] = min_value  # Replace P=0 with min non-0 p-value
    enrich_ps['-log10(p-value)'] = abs(-1 * np.log10(enrich_ps[pval_col]))

    enrich_ps.reset_index(drop=True, inplace=True)

    return enrich_ps


def plot_interactive_phewas(data, x_column='trait_reported',
                            y_column='-log10(p-value)',
                            color_column='trait_category',
                            filter_column='program_name',
                            significance_threshold=0.05,
                            annotation_cols=["program_name", "trait_reported",
                                             "trait_category", "P-value",
                                             "Genes", "study_id", "pmid"],
                           query_string="trait_category != 'measurement'",
                           title="Cell Program x OpenTargets GWAS L2G Enrichment"):
    
    # Get unique values for the filtering column
    filter_values = ['All'] + list(data[filter_column].unique())

    # Initialize output widget to display the plot
    output = Output()
    
    if query_string:
        data=data.query(query_string)

    # Function to update plot based on dropdown selection
    def update_plot(selected_value):
        # Filter data based on selected value
        if selected_value == "All":
            filtered_data = data.copy()  # No selection, show all data
        else:
            filtered_data = data[data[filter_column] == selected_value]

        # Create the plot
        fig = px.scatter(filtered_data, x=x_column, y=y_column, color=color_column,
                         title=title,
                         hover_data=annotation_cols)

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