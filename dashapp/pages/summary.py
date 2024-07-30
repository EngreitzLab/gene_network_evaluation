import dash
from dash import html, dcc, callback, Input, Output


def summary_page():
    return html.Div(
        id='summary_page', 
        className='section', 
        style={'width': '50%', 'display': 'inline-block', 'vertical-align': 'top'}, 
        children=[
            html.H2("Input Data Summary"),
            html.P("Number of AnnData’s analyzed: ..."),  # Placeholder
            html.P("Number of cells in each AnnData: ..."),  # Placeholder
            html.P("This section will be filled with actual data summary."),
        ]
)


dash.register_page(__name__)

layout = summary_page()
