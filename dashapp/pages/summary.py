import dash
from dash import html, dcc, callback, Input, Output


def summary_page():
    return html.Div(id='landing_page', className='section', style={'width': '50%', 'display': 'inline-block', 'vertical-align': 'top'}, children=[
        html.H2("Input Data Summary"),
        html.P(f"{len(mdata.mod)} AnnData’s analyzed: {', '.join(mdata.mod.keys())}"),
        html.P(f"Number of cells in each AnnData: {summary_dict['n_cells']}"),
        html.P("Number of programs analyzed per modality:"),
        html.Ul([html.Li(f"{mod}: {summary_dict[mod]['n_programs']}") for mod in mdata.mod.keys()]),
    ])


dash.register_page(__name__)

layout = summary_page()
