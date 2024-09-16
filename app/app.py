import os
import sys
import yaml
import logging
import argparse
import collections
import mudata
import pandas as pd
import dash
from dash import Dash, html, dcc
import dash_bootstrap_components as dbc
import diskcache
from dash import DiskcacheManager
from uuid import uuid4
from parse import parse
from utils import infer_dashboard_type, load_config
    

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
        path_evaluation_outs,
        data_key="rna",
        perturbation_association_stratification_key=None,
        motif_enrichment_stratification_key=None,
    ):
        try:
            # Load mdata
            mdata = mudata.read_h5mu(path_mdata)
            mdata.mod = collections.OrderedDict(sorted(mdata.mod.items()))
            
            # Parse data
            results = parse(
                mdata, 
                path_evaluation_outs, 
                data_key,
                perturbation_association_stratification_key=perturbation_association_stratification_key,
                motif_enrichment_stratification_key=motif_enrichment_stratification_key,
            )
            
            # Add obs data from data_key
            obs = mdata[data_key].obs.reset_index()
            obs.columns = ["barcode"] + list(obs.columns[1:])
            results["obs"] = obs
            
            # Add obsm data
            obsms = {}
            for obsm_key in mdata[data_key].obsm:
                cols = [f"{obsm_key}_{i}" for i in range(2)]
                rows = mdata[data_key].obs_names
                df = pd.DataFrame(mdata[data_key].obsm[obsm_key][:, :2], columns=cols, index=rows)
                obsms[obsm_key] = df
            results["obsms"] = obsms

            # Return results
            return results

        except Exception as e:
            print(f"Could not load mdata: {e}")
            sys.exit(1)

        
    # Parse config for paths
    path_evaluation_outs = config["path_evaluation_outs"]
    path_mdata = config["path_mdata"]
    path_evaluation_config = config["path_evaluation_config"]
    path_report_out = config["path_report_out"]

    # Parse config for other parameters
    data_key = config["data_key"]
    categorical_keys = config["categorical_keys"] if config["categorical_keys"] else []
    continuous_keys = config["continuous_keys"] if config["continuous_keys"] else []
    annotations_loc = config["annotations_loc"]

    # Load evaluation config
    evaluation_config = load_config(path_evaluation_config)
    perturbation_association_stratification_key = evaluation_config["perturbation_association"]["groupby_key"]
    motif_enrichment_stratification_key = evaluation_config["motif_enrichment"]["groupby_key"]
    
    # Set up logging to print to console and also to file in path_out (evaluation_pipeline.log) with overwrite
    log_path = os.path.join(path_report_out, 'reports_pipeline.log')
    if os.path.exists(log_path):
        os.remove(log_path)
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', handlers=[
        logging.FileHandler(log_path, mode='w'),
        logging.StreamHandler()
    ])

    # Log configuration
    for key, value in config.items():
        logging.info(f"{key}: {value}")

    # Load and parse data
    results = load_and_parse_data(
        path_mdata=path_mdata,
        path_evaluation_outs=path_evaluation_outs,
        data_key=data_key,
        perturbation_association_stratification_key=perturbation_association_stratification_key,
        motif_enrichment_stratification_key=motif_enrichment_stratification_key,
    )
    
    # Add in path_report_out to results
    results['path_evaluation_outs'] = path_evaluation_outs
    results['path_mdata'] = path_mdata
    results['evaluation_config'] = evaluation_config
    results['path_report_out'] = path_report_out
    results['data_key'] = data_key
    results['categorical_keys'] = categorical_keys
    results['continuous_keys'] = continuous_keys
    results['annotations_loc'] = annotations_loc
    results['perturbation_association_stratification_key'] = perturbation_association_stratification_key
    results['motif_enrichment_stratification_key'] = motif_enrichment_stratification_key

    # Cache results
    cache.set('results', results)
    print(f"Loaded main data: {results.keys()}")
    
    # Get program keys from methods in results
    prog_keys = list(results["methods"].keys())
    print(f"Program keys: {prog_keys}")

    # Dashboard 
    dashboard_type = infer_dashboard_type(prog_keys)
    print(f"Dashboard type: {dashboard_type}")

    # Create Dash app
    app = Dash(
        __name__, 
        use_pages=True, 
        background_callback_manager=background_callback_manager, 
        external_stylesheets=[dbc.themes.LUX],  # Updated theme to LUX for a modern look
        suppress_callback_exceptions=True
    )
    app.title = "Gene Program Evaluation Dashboard v0.0.1"

    # Define Sidebar
    sidebar = dbc.Nav(
        [
            dbc.NavLink(
                [
                    html.Div(page["name"], className="sidebar-link"),
                ],
                href=page["relative_path"],
                active="exact"
            )
            for page in dash.page_registry.values()
        ],
        vertical=True,
        pills=True,
        className="sidebar",
    )

    # Define the app layout with modern, minimalistic style
    app.layout = dbc.Container([
        
        # Title Row
        dbc.Row([
            dbc.Col(html.H1("Gene Program Evaluation Dashboard", 
                            style={'fontSize': '3rem', 'textAlign': 'center', 
                                   'color': '#2c3e50', 'marginBottom': '20px', 'fontWeight': 'bold'}))
        ]),

        html.Hr(style={'borderColor': '#2c3e50'}),  # Horizontal rule with a nice color

        # Main Layout with Sidebar and Content
        dbc.Row(
            [
                # Sidebar
                dbc.Col(
                    [
                        sidebar
                    ], xs=4, sm=4, md=3, lg=2, xl=2, className="bg-light sidebar-container"
                ),

                # Page Content
                dbc.Col(
                    [
                        dash.page_container
                    ], xs=8, sm=8, md=9, lg=10, xl=10
                )
            ],
            className="main-content"
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
    path_report_out = config["path_report_out"]
    os.makedirs(path_report_out, exist_ok=True)
    app = create_dash_app(config)
    app.run_server(debug=True, host='0.0.0.0', port=8050)
