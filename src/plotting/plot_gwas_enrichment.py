import pandas as pd
import numpy as np
import mudata
import plotly.express as px
from ipywidgets import interact, Dropdown, Output, VBox
from IPython.display import display

import pandas as pd
import numpy as np


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