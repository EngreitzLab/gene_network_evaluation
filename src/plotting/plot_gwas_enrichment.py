import pandas as pd
import numpy as np
import mudata
import plotly.express as px
from ipywidgets import interact, Dropdown, Output, VBox

def plot_interactive_phewas(data, x_column='trait_reported', y_column='-log10(p-value)', color_column='trait_category', filter_column='Program_Name'):
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
            height=500,  # Adjust height as needed
            xaxis_tickfont=dict(size=4)
        )

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