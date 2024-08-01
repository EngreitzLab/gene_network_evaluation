import dash
from dash import Dash, html, dcc, Input, Output
from data_processing import parse_mdata_summary, parse_mdata_cross_run, extract_total_unique_counts
from layout import create_scatter_layout, create_filtered_barplot_layout, create_filtered_stacked_barplot_layout
import dash_bootstrap_components as dbc
import mudata


# Load mudata
mdata = mudata.read_h5mu("./example_data/Endothelial/cNMF_evaluation_dashapp_data_small.h5mu")
summary_dict = parse_mdata_summary(mdata, verbose=False)
cross_run_dict = parse_mdata_cross_run(mdata, verbose=False)

# Create Dash app
app = Dash(__name__, use_pages=True, pages_folder="./pages", external_stylesheets=[dbc.themes.SPACELAB])
app.title = "Gene Program Evaluation Dashboard v0.0.1"

app.layout = html.Div(
    [
        dcc.Store(id="store", data={}),
        html.H1("Gene Program Evaluation Dashboard v0.0.1"),
        html.Div(
            [
                html.Div(
                    dcc.Link(f"{page['name']}", href=page["path"]),
                )
                for page in dash.page_registry.values()
            ]
        ),
        html.Hr(),
        dash.page_container,
    ]
)


# Run the app
if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0', port=8050)
