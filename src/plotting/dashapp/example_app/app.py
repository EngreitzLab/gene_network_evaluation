from dash import Dash
from data_processing import load_data, count_unique
from layout import create_layout
from callbacks import register_callbacks

# Load data
explained_variance, enrichment_gsea, enrichment_motif, enrichment_trait = load_data('/cellar/users/aklie/opt/gene_program_evaluation/src/plotting/dashapp/example_data/cNMF_evaluation_output.xlsx')

# Preprocess data
gsea_unique_df = count_unique(
    categorical_var='ProgramID',
    count_var='ID',
    dataframe=enrichment_gsea.loc[enrichment_gsea['qvalue'] <= 0.05]
)

# Create Dash app
app = Dash(__name__)
app.layout = create_layout(explained_variance, gsea_unique_df, enrichment_gsea, enrichment_motif, enrichment_trait)

# Register callbacks
register_callbacks(app, enrichment_gsea, enrichment_motif, enrichment_trait)

if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0', port=8050)
