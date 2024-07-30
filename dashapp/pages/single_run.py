import dash
from dash import html, dcc, callback, Input, Output

def single_run_analysis():
    return html.Div(id='single_run_analysis', className='section', style={'width': '50%', 'display': 'inline-block', 'vertical-align': 'top'}, children=[
        html.H2("Single Run Analysis"),
        html.P("Description: Analyze individual runs."),
        html.P("This section will contain detailed analysis for a selected run.")
    ])

dash.register_page(__name__)

layout = single_run_analysis()
