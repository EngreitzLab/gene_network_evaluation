import pandas as pd
import muon as mu
import scanpy as sc
import celloracle as co
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-d','--path_data', required=True)
parser.add_argument('-r','--path_r2g', required=True)
parser.add_argument('-t','--path_tf2r', required=True)
parser.add_argument('-c','--cluster_key', required=True)
parser.add_argument('-l','--layer', required=False)
parser.add_argument('-a','--alpha', required=True)
parser.add_argument('-b','--bagging_number', required=True)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

# Load args
path_data = args['path_data']
path_r2g = args['path_r2g']
path_tf2r = args['path_tf2r']
cluster_key = args['cluster_key']
layer = args['layer']
alpha = float(args['alpha'])
bagging_number = int(args['bagging_number'])
path_out = args['path_out']

# Process base GRN
r2g = pd.read_csv(path_r2g)
tfb = pd.read_csv(path_tf2r)
if (r2g.shape[0] == 0) or (tfb.shape[0] == 0):
    grn = pd.DataFrame(columns=['source', 'target', 'score', 'pval'])
    grn.to_csv(path_out, index=False)
    exit()
tfb['score'] = 1
r2g = r2g[['cre', 'gene']]
base_grn = pd.merge(
    r2g,
    tfb
    .pivot(index='cre', columns='tf')
    .fillna(0)
    .droplevel(0, axis=1)
    .reset_index()
)
base_grn = base_grn.rename(columns={'cre': 'peak_id', 'gene': 'gene_short_name'})
base_grn['peak_id'] = base_grn['peak_id'].str.replace('-', '_')

# Init object
mdata = mu.read(path_data)
adata = mdata.mod["rna"].copy()
adata.obs[cluster_key] = mdata.obs[cluster_key].copy()
if layer is not None:
    print(f"Using data in layer {layer} for regression.")
    adata.X = adata.layers[layer].copy()
else:
    print("Log normalizing data for regression.")
    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

# Model TF ~ G for every cluster
cluster_grns = {}
for cluster in adata.obs[cluster_key].cat.categories:
    print(f"Building GRN for {cluster}")
    adata_sub = adata[adata.obs[cluster_key] == cluster].copy()
    net = co.Net(
        gene_expression_matrix=adata_sub.to_df(), # Input gene expression matrix as data frame
        TFinfo_matrix=base_grn, # Input base GRN
        verbose=True
    )
    net.fit_All_genes(
        bagging_number=bagging_number,
        alpha=alpha,
        verbose=True
    )
    net.updateLinkList(verbose=True)
    inference_result = net.linkList.copy()
    cluster_grns[cluster] = inference_result
    print(f"Finished building GRN for {cluster}")

# Extract grn
grn = pd.concat([v.assign(cluster=k) for k, v in cluster_grns.items()])
grn = grn.dropna()[['source', 'target', 'coef_mean', 'p', 'cluster']]
grn = grn.rename(columns={'coef_mean': 'score', 'p': 'pval'})
grn = grn.sort_values(['source', 'score'], ascending=[True, False])

# Write
grn.to_csv(path_out, index=False)
