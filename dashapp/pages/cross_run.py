import os
import dash
import pickle
from dash import html, dcc, callback, Input, Output

# Ouput directory
path_pipeline_outs = "/cellar/users/aklie/opt/gene_program_evaluation/dashapp/example_data/iPSC_EC_evaluations"
data_key = "rna"

dash.register_page(__name__, order=1)

# Load results from pickle
with open(os.path.join(path_pipeline_outs, "results.pkl"), "rb") as f:
    results = pickle.load(f)

print(results.keys())

layout = html.Div(id='cross_run', className='section', style={'width': '50%', 'display': 'inline-block', 'vertical-align': 'top'}, children=[
    html.H2("Cross Run"),
    html.Div("This is the cross run page."),
])