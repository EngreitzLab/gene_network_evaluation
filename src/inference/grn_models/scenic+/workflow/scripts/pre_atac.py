import os
import pickle
import mudata as mu
from pycisTopic.cistopic_class import create_cistopic_object
from pycisTopic.lda_models import run_cgs_models, evaluate_models
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.diff_features import impute_accessibility, normalize_scores, find_highly_variable_features, find_diff_features
from pycisTopic.utils import region_names_to_coordinates
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_data', required=True)
parser.add_argument('-p','--path_out', required=True)
parser.add_argument('-c','--cluster_key', required=True)
parser.add_argument('-o','--organism', required=True)
parser.add_argument('-n','--n_topics', required=True)
parser.add_argument('-t','--n_iter', required=True)
parser.add_argument('-a','--alpha', required=True)
parser.add_argument('-x','--alpha_by_topic', required=True)
parser.add_argument('-e','--eta', required=True)
parser.add_argument('-y','--eta_by_topic', required=True)
parser.add_argument('-nt','--ntop', required=True)
parser.add_argument('-si','--scale_factor_impute', required=True)
parser.add_argument('-sn','--scale_factor_normalize', required=True)
parser.add_argument('-j','--n_cpu', required=True)
parser.add_argument('-d','--temp_dir', required=True)
parser.add_argument('-r','--random_state', required=True)
args = vars(parser.parse_args())

# Parse args
path_data = args['path_data']
path_out = args['path_out']
cluster_key = args['cluster_key']
organism = args['organism']
n_topics = args['n_topics']
n_iter = int(args['n_iter'])
alpha = int(args['alpha'])
alpha_by_topic = bool(args['alpha_by_topic'])
eta = float(args['eta'])
eta_by_topic = bool(args['eta_by_topic'])
ntop = int(args['ntop'])
scale_factor_impute = float(args['scale_factor_impute'])
scale_factor_normalize = float(args['scale_factor_normalize'])
n_cpu = int(args['n_cpu'])
temp_dir = args['temp_dir']
random_state = int(args['random_state'])

# Read atac adata
adata = mu.read(path_data)
del adata.mod['rna']
obs = adata.obs.copy()
adata = adata.mod['atac'].copy()
adata.obs = obs

# Clean up cluster_key if necessary 
adata.obs[cluster_key] = adata.obs[cluster_key].str.replace("/", "_").copy()

if organism == 'human':
    path_blacklist = 'resources/blacklists/human.bed'
elif organism == 'mouse':
    path_blacklist = 'resources/blacklists/mouse.bed'

# Make first occurences of "-" ":" in var_names if multiple exist
if sum(adata.var.index.str.contains("-")) > 0:
    adata.var.index = adata.var.index.str.replace("-", ":", 1)

# Create cisTopic object
cistopic_obj = create_cistopic_object(
    fragment_matrix=adata.to_df().T,
    cell_names=adata.obs.index.values,
    region_names=adata.var.index.values,
    path_to_blacklist=path_blacklist,
    split_pattern="_",
    tag_cells=False,
)

# Add cell metadata
cistopic_obj.add_cell_data(obs)

# Process n_topics
n_topics = n_topics.split(';')
n_topics = [int(n) for n in n_topics]

# Compute topic models
models = run_cgs_models(
    cistopic_obj,
    n_topics=n_topics,
    n_cpu=n_cpu,
    n_iter=n_iter,
    random_state=random_state,
    alpha=alpha,
    alpha_by_topic=alpha_by_topic,
    eta=eta,
    eta_by_topic=eta_by_topic,
    save_path=None,
    _temp_dir=temp_dir
)

# Save models object
with open(os.path.join(path_out, "models.pkl"), 'wb') as f:
    pickle.dump(models, f)

# Add model to cistopic object
model = evaluate_models(
    models,
    select_model=None,
    return_model=True,
    metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
    plot_metrics=False
)
cistopic_obj.add_LDA_model(model)

# Save the cistopic object
pickle.dump(cistopic_obj, open(os.path.join(path_out, "cistopic_obj.pkl"), 'wb'))

# Binarize topics using otsu and top
region_bin_topics_otsu = binarize_topics(
    cistopic_obj, 
    method='otsu',
    plot=False
)
region_bin_topics_top = binarize_topics(
    cistopic_obj, 
    method='ntop', 
    ntop=ntop,
    plot=False
)

# Impute accessibility by matrix multiplication
imputed_acc_obj = impute_accessibility(
    cistopic_obj,  # cisTopic object
    selected_cells=None,  # certain cells to use for imputation
    selected_regions=None,  # certain regions to use for imputation
    scale_factor=scale_factor_impute,  # A number to multiply the imputed values for.
    project=cistopic_obj.project + "impute"  # A string to add to the name of the imputed object.
)

# Normalize accessibility scores
normalized_imputed_acc_obj = normalize_scores(
    imputed_acc_obj, 
    scale_factor=scale_factor_normalize  # Scale factor for cell-level normalization
)

# Find highly variable regions
variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj, 
    plot=False
)
del normalized_imputed_acc_obj

# Find differentially accessible regions
markers = find_diff_features(
    cistopic_obj, 
    imputed_acc_obj, 
    variable=cluster_key,
    var_features=variable_regions, 
    split_pattern="_"
)

# Make dirs for saving region sets
os.makedirs(os.path.join(path_out, "region_sets"), exist_ok = True)
os.makedirs(os.path.join(path_out, "region_sets", "Topics_otsu"), exist_ok = True)
os.makedirs(os.path.join(path_out, "region_sets", "Topics_top"), exist_ok = True)
os.makedirs(os.path.join(path_out, "region_sets", "DARs"), exist_ok = True)

# Save region sets (TODO: Make function)
for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(path_out, "region_sets", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )
for topic in region_bin_topics_top:
    region_names_to_coordinates(
        region_bin_topics_top[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(path_out, "region_sets", "Topics_top", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )
for cluster in markers:
    region_names_to_coordinates(
        markers[cluster].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(path_out, "region_sets", "DARs", f"{cluster}.bed"),
        sep = "\t",
        header = False, index = False
    )
