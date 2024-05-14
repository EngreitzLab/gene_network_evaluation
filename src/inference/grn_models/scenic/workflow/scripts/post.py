import os
import glob
import numpy as np
import mudata as mu
import scanpy as sc
import loompy as lp
import pandas as pd
from pyscenic.cli.utils import load_signatures
from scipy.stats import ttest_1samp
from pyscenic.utils import add_correlation
from tqdm import tqdm
import argparse

# 
tqdm.pandas()
TINY = np.finfo(np.float32).tiny

# Init args
parser = argparse.ArgumentParser(
    prog="python post.py",
    description="Post process outputs of pipeline for evaluation."
)
parser.add_argument('-i','--path_data', required=True)
parser.add_argument('-l','--path_loom', required=True)
parser.add_argument('-c','--path_csvs', required=True, nargs='+')
parser.add_argument('-o','--path_output', required=True)
args = vars(parser.parse_args())

# Parse args
path_data = args['path_data']
path_loom = args['path_loom']
path_csvs = args['path_csvs']
path_output = args['path_output']

# functions
def calc_p_value(importances):
    _, p_value = ttest_1samp(importances, 0)
    return p_value

def regulon2sadj(
    regulons,
):
    net_lst = []
    for tf in regulons:
        tf_name = tf.name.split("(")[0]
        tf_targets = tf.gene2weight
        for target, weight in tf_targets.items():
            net_lst.append([tf_name, target, weight])
    net = pd.DataFrame(net_lst, columns=["TF", "target", "importance"])
    return net

# Read mdata to add objects to
data = mu.read(path_data)

# Read regulons into df
print("Reading regulons...")
all_edges = pd.DataFrame()
for reg_csv in path_csvs:
    regulons = load_signatures(reg_csv)
    adj_df = regulon2sadj(regulons)
    all_edges = pd.concat([all_edges, adj_df])
all_edges.head()
print(f"Total edges: {len(all_edges)}")


# Post process
print("Grouping by source and target and filtering singlet edges...")
grouped = all_edges.groupby(['TF', 'target'])
filtered = grouped.filter(lambda x: len(x) > 1)
print(f"{len(all_edges) - len(filtered)} edges dropped")

print("Calculating mean importance for each edge...")
mean_importance = filtered.groupby(['TF', 'target'])['importance'].mean()
print(f"Total unique edges: {len(mean_importance)}")

print("Calculating empirical p-value...")
p_values_series = filtered.groupby(['TF', 'target'])['importance'].progress_apply(calc_p_value)
p_values = p_values_series.values + TINY

print("Transforming values...")
neg_log_p = -np.log10(p_values)
normalized_importance = (mean_importance - mean_importance.min()) / (mean_importance.max() - mean_importance.min())

print("Adding correlation...")
adata = sc.read_loom(path_loom, sparse=True)
filtered = add_correlation(filtered, adata.to_df())
mean_corr = filtered.groupby(['TF', 'target'])['rho'].mean()

# Consolidate
consolidated = pd.DataFrame({
    'tf': mean_importance.index.get_level_values('TF'),
    'gene': mean_importance.index.get_level_values('target'),
    'weight_signed': np.nan,
    'weight_unsigned': mean_importance.values,
    'weight_minmax_normalized': normalized_importance.values,
    'pval': p_values,
    '-logp': neg_log_p,
    'description': np.nan,
    'corr': mean_corr.values,
    'cluster': 'global'
}).reset_index(drop=True)
consolidated.to_csv(os.path.join(path_csvs, "full_grn.csv"))

# Remove self-loops
print("Removing self-loops...")
consolidated = consolidated[consolidated["tf"] != consolidated["gene"]]

# GRN
grn = consolidated[["tf", "gene", "weight_unsigned", "pval", "cluster"]]
grn = grn.rename({"weight_unsigned": "score"}, axis=1)
grn.to_csv(os.path.join(path_csvs, "grn.csv"))

# add to data
data.uns["grn"] = grn

# Save data
data.write_h5mu(path_output)
