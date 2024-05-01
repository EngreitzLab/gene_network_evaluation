import numpy as np
import scanpy as sc
import celloracle as co
import muon as mu
import os
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_input', required=True)
parser.add_argument('-k','--knn', required=True)
parser.add_argument('-d','--dim_reduction_key', required=False)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

path_input = args['path_input']
k = int(args['knn'])
dim_reduction_key = args['dim_reduction_key']
path_out = args['path_out']

# Read rna adata
mdata = mu.read(path_input)

# Extract raw counts data and assign labels
adata = mdata.mod['rna'].copy()
adata.obs['celltype'] = mdata.obs['celltype']
adata.X = adata.layers['counts'].copy()

# Get dim reduction if none passed in
if dim_reduction_key is None:
    print("No dim reduction key passed in. Using PCA.")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=50)
    adata.X = adata.layers['counts'].copy()
    dim_reduction_key = "X_pca"

# Instantiate Oracle object
oracle = co.Oracle()
oracle.import_anndata_as_raw_count(
    adata=adata,
    cluster_column_name="celltype",
    embedding_name=dim_reduction_key,
)

# Compute PCA and select top pcs
oracle.perform_PCA()
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
n_comps = min(n_comps, 50)

# Run imputation
oracle.knn_imputation(
    n_pca_dims=n_comps,
    k=k,
    balanced=True,
    b_sight=k*8,
    b_maxl=k*4,
    n_jobs=os.cpu_count(),
)

# Update object with imputet counts
mdata['rna'].X = oracle.adata.layers['imputed_count']

# Write
mdata.write(path_out)
