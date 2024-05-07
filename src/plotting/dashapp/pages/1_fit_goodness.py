import numpy as np
import pandas as pd

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from dash import Dash, html, dcc

import argparse

# Plot goodness of fit heuristics
def run_app(fit_metrics: pd.DataFrame, num_cols=3):

    app = Dash(__name__)

    num_metrics = len(fit_metrics.columns) - 1
    num_rows = int(num_metrics/num_cols) + int(np.ceil(num_metrics%num_cols))

    prog_name_col = [col for col in fit_metrics.columns if 'prog_name' in col][0]

    fig = make_subplots(
        rows=num_rows, cols=num_cols,
        specs=[[{}]*num_cols]*num_rows,
        print_grid=True,
        subplot_titles=['{}'.format(col) for col in \
                        fit_metrics.columns if 'prog_name' not in col],
        vertical_spacing = 0.05, horizontal_spacing = 0.1)

    plots = {}
    for i, col in enumerate([col for col in fit_metrics.columns if 'prog_name' not in col]):

        plots[col] = px.scatter(x=fit_metrics[prog_name_col], 
                                y=fit_metrics[col])
        plots[col].update_layout(xaxis_title='Components', yaxis_title=col)
        
        row_num = int((i)/num_cols) + 1
        col_num = i - (row_num-1)*num_cols + 1
        fig.add_trace(plots[col]['data'][0], row=row_num, col=col_num)
        fig.update_xaxes(showticklabels=False, row=row_num, col=col_num)
        fig.update_yaxes(ticksuffix = "  ", row=row_num, col=col_num)

    fig.update_traces(hovertemplate="%{y}")
    fig.update_layout(hovermode="x unified")

    app.layout = html.Div(children=[
        html.H1(children='GEP Dashboard - v0.1"'),

        html.Div(children='''
            Goodness of fit measures to infer optimal number of programs.
              '''),

        dcc.Graph(
            id='gof_measures',
            figure=fig
        )
    ])

    app.run(debug=True)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('fit_metrics_path', type=str)

    args = parser.parse_args()
    
    fit_metrics = pd.read_csv(args.fit_metrics_path, sep='\t')
    run_app(fit_metrics)
