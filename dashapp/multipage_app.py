import os
import sys
import pickle
import dash
from dash import Dash, html, dcc, Input, Output
import mudata
import dash_bootstrap_components as dbc
import pandas as pd
import collections
from parse import parse
from utils import infer_dashboard_type

# Ouput directory
path_pipeline_outs = "/cellar/users/aklie/opt/gene_program_evaluation/dashapp/example_data/iPSC_EC_evaluations"
data_key = "rna"

# Load mdata
try:
    path_mdata = os.path.join(path_pipeline_outs, "cNMF_60_0.2_gene_names.h5mu")
    mdata = mudata.read_h5mu(path_mdata)
    mdata.mod = collections.OrderedDict(sorted(mdata.mod.items()))
    anndata_keys = list(mdata.mod.keys())
    prog_keys = [key for key in anndata_keys if key != data_key]
    with mudata.set_options(display_style="html", display_html_expand=0b000):
        html_rep = mdata._repr_html_(expand=0b000)
    with open(os.path.join(path_pipeline_outs, "mudata.html"), "w") as f:
        f.write(html_rep)
    for obsm_key in mdata[data_key].obsm:
        cols = [f"{obsm_key}_{i}" for i in range(2)]
        rows = mdata[data_key].obs_names
        df = pd.DataFrame(mdata[data_key].obsm[obsm_key][:, :2], columns=cols, index=rows)
        # Add columns for all .obs categorical columns
        for col in mdata.obs.columns:
            if pd.api.types.is_categorical_dtype(mdata.obs[col]):
                df[col] = mdata.obs[col].values
        df.to_csv(os.path.join(path_pipeline_outs, f"{obsm_key}.tsv"), sep="\t", index=True)
except:
    print("Could not load mdata")
    sys.exit(1)

# Get all subdirectories of the pipeline outputs
subdirs = [x[0] for x in os.walk(path_pipeline_outs)][1:]

# Parse and save as pickle
results = parse(mdata, subdirs, data_key)
with open(os.path.join(path_pipeline_outs, "results.pkl"), "wb") as f:
    pickle.dump(results, f)

# Dashboard 
dashboard_type = infer_dashboard_type(prog_keys)
print(f"Dashboard type: {dashboard_type}")

# Create Dash app
app = Dash(__name__, use_pages=True, external_stylesheets=[dbc.themes.SPACELAB], suppress_callback_exceptions=True)
app.title = "Gene Program Evaluation Dashboard v0.0.1"

# If the dashboard type is single_run, then we want to remove pages.cross_run
if dashboard_type == "single_run":
    del dash.page_registry["pages.cross_run"]

# Create sidebar for page navigation
sidebar = dbc.Nav(
    [
        dbc.NavLink(
            [
                html.Div(page["name"], className="sidebar-link"),
            ],
            href=page["relative_path"],
        )
        for page in dash.page_registry.values()
    ],
    vertical=True,
    pills=True,
    className="sidebar",
)

app.layout = dbc.Container([
    dbc.Row([
        dbc.Col(html.Div("Gene Program Evaluation Dashboard v0.0.1",
                         style={'fontSize':50, 'textAlign':'center'}))
    ]),

    html.Hr(),

    dbc.Row(
        [
            dbc.Col(
                [
                    sidebar
                ], xs=4, sm=4, md=2, lg=2, xl=2, xxl=2),

            dbc.Col(
                [
                    dash.page_container
                ], xs=8, sm=8, md=10, lg=10, xl=10, xxl=10)
        ]
    )
], fluid=True)

# Run the app
if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0', port=8050)
