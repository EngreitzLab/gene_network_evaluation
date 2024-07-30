import dash
from dash import html, dcc, callback, Input, Output

def program_analysis():
    return html.Div(id=f'program_analysis', className='section', style={'width': '50%', 'display': 'inline-block', 'vertical-align': 'top'}, children=[
        html.H2("Program Analysis"),
        html.P("Description: Analyze individual programs."),
        html.P("This section will contain detailed analysis for a selected program.")
    ])


dash.register_page(__name__)

layout = program_analysis()
