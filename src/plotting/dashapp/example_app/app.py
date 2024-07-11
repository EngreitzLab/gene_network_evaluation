from dash import Dash
from data_processing import load_data, count_unique
from layout import create_layout
import pandas as pd
from callbacks import register_callbacks

# Load data
explained_variance, enrichment_gsea, enrichment_motif, enrichment_trait = load_data('../example_data/cNMF_evaluation_output.xlsx')
trait_phewas = pd.read_csv("../example_data/cNMF_enrichment_trait_processed.txt", sep='\t')

# Create Dash app
app = Dash(__name__)
app.layout = create_layout(explained_variance, enrichment_gsea, enrichment_motif, enrichment_trait, trait_phewas)

# Register callbacks
register_callbacks(app, enrichment_gsea, enrichment_motif, enrichment_trait, trait_phewas)

if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0', port=8050)
