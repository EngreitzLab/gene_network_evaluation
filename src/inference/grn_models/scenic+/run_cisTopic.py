import pickle
import pandas as pd
import numpy as np
import muon as mu
import argparse
from pycisTopic.cistopic_class import create_cistopic_object
from pycisTopic.lda_models import run_cgs_models, evaluate_models


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_input', required=True)
parser.add_argument('-o','--organism', required=True)
parser.add_argument('-n','--n_topics', required=True)
parser.add_argument('-t','--n_iter', required=True)
parser.add_argument('-a','--alpha', required=True)
parser.add_argument('-x','--alpha_by_topic', required=True)
parser.add_argument('-e','--eta', required=True)
parser.add_argument('-y','--eta_by_topic', required=True)
parser.add_argument('-m','--path_models', required=True)
parser.add_argument('-c','--path_cistopic_obj', required=True)
args = vars(parser.parse_args())

path_input = args['path_input']
organism = args['organism']
path_cistopic_obj = args['path_cistopic_obj']
n_topics = args['n_topics']
n_iter = int(args['n_iter'])
alpha = int(args['alpha'])
alpha_by_topic = bool(args['alpha_by_topic'])
eta = float(args['eta'])
eta_by_topic = bool(args['eta_by_topic'])
path_models = args['path_models']

# Read atac adata
adata = mu.read(path_input)
del adata.mod['rna']
obs = adata.obs.copy()
adata = adata.mod['atac'].copy()
adata.obs = obs

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
    n_cpu=1,  # TODO: no hardcoding
    n_iter=n_iter,
    random_state=2017,
    alpha=alpha,
    alpha_by_topic=alpha_by_topic,
    eta=eta,
    eta_by_topic=eta_by_topic,
    save_path=None,
    _temp_dir="/cellar/users/aklie/tmp",  # TODO: no hardcoding
)

# Save models object
with open(path_models, 'wb') as f:
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
pickle.dump(cistopic_obj, open(path_cistopic_obj, 'wb'))
