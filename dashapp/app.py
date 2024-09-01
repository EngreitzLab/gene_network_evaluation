import os
import sys
import yaml
import logging
import argparse
import collections
import mudata
import pandas as pd
import dash
from dash import Dash, html
import dash_bootstrap_components as dbc
import diskcache
from dash import DiskcacheManager
from uuid import uuid4
from parse import parse
from utils import infer_dashboard_type

def load_config(config_path):
    with open(config_path, 'r') as config_file:
        return yaml.safe_load(config_file)
    
def create_dash_app(config):
    
    # Setup cache
    launch_uid = uuid4()
    cache = diskcache.Cache("./.cache")
    background_callback_manager = DiskcacheManager(
        cache, cache_by=[lambda: launch_uid], expire=60
    )

    # Cache the data loading and parsing function
    @cache.memoize()
    def load_and_parse_data(
        path_mdata,
        path_pipeline_outs,
        data_key="rna",
        categorical_keys=[],
    ):
        try:
            
            # Load mdata
            mdata = mudata.read_h5mu(path_mdata)
            mdata.mod = collections.OrderedDict(sorted(mdata.mod.items()))

            # Get subdirectories
            subdirs = [x[0] for x in os.walk(path_pipeline_outs)][1:]
            
            # Parse data
            results = parse(mdata, subdirs, data_key)
            
            # Add MuData HTML representation
            with mudata.set_options(display_style="html", display_html_expand=0b000):
                html_rep = mdata._repr_html_(expand=0b000)
            results["mdata"] = html_rep
            
            # Add obsm data
            obsms = {}
            for obsm_key in mdata[data_key].obsm:
                cols = [f"{obsm_key}_{i}" for i in range(2)]
                rows = mdata[data_key].obs_names
                df = pd.DataFrame(mdata[data_key].obsm[obsm_key][:, :2], columns=cols, index=rows)
                for col in mdata.mod[data_key].obs.columns:
                    if col in categorical_keys:
                        if pd.api.types.is_categorical_dtype(mdata[data_key].obs[col]):
                            df[col] = mdata[data_key].obs[col].values
                obsms[obsm_key] = df
            results["obsms"] = obsms

            # Return results
            return results

        except Exception as e:
            print(f"Could not load mdata: {e}")
            sys.exit(1)

        
    # Parse config for paths
    path_out = config["path_out"]
    path_mdata = config["path_mdata"]
    data_key = config["data_key"]
    categorical_keys = config["categorical_keys"]
    os.makedirs(path_out, exist_ok=True)
    
    # Set up logging to print to console and also to file in path_out (evaluation_pipeline.log) with overwrite
    log_path = os.path.join(path_out, 'reports_pipeline.log')
    if os.path.exists(log_path):
        os.remove(log_path)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', handlers=[
        logging.FileHandler(log_path, mode='w'),
        logging.StreamHandler()
    ])

    # Log configuration
    logging.info(f"Configuration: {config}")

    # Load and parse data
    results = load_and_parse_data(
        path_mdata=path_mdata,
        path_pipeline_outs=path_out,
        data_key=data_key,
        categorical_keys=categorical_keys
    )
    cache.set('results', results)
    print(f"Loaded main data: {results.keys()}")
    
    # Get program keys from methods in results
    prog_keys = list(results["methods"].keys())
    print(f"Program keys: {prog_keys}")

    # Get all subdirectories of the pipeline outputs
    subdirs = [x[0] for x in os.walk(path_out)][1:]
    print(f"Subdirectories: {subdirs}")

    # Dashboard 
    dashboard_type = infer_dashboard_type(prog_keys)
    print(f"Dashboard type: {dashboard_type}")

    # Create Dash app
    app = Dash(
        __name__, 
        use_pages=True, 
        background_callback_manager=background_callback_manager, 
        external_stylesheets=[dbc.themes.SPACELAB], 
        suppress_callback_exceptions=True
    )
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

    # Define app layout including dcc.Store for sharing data across pages
    app.layout = dbc.Container([
        
        dbc.Row([
            dbc.Col(html.Div("Gene Program Evaluation Dashboard v0.0.1",
                            style={'fontSize': 50, 'textAlign': 'center'}))
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

    return app

# Run the app
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run the Dash app with a configuration file.")
    parser.add_argument(
        '--config',
        type=str,
        required=True,
        help="Path to the configuration file."
    )
    args = parser.parse_args()
    config = load_config(args.config)
    app = create_dash_app(config)
    app.run_server(debug=True, host='0.0.0.0', port=8050)
