import os
import dash
import pickle
from dash import html, dcc, callback, Input, Output


dash.register_page(__name__, order=1)

layout = html.Div(id='cross_run', className='section', style={'width': '50%', 'display': 'inline-block', 'vertical-align': 'top'}, children=[
    html.H2("Cross Run"),
    html.Div("This is the cross run page."),
])