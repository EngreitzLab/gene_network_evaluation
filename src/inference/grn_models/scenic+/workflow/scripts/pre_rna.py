import mudata as mu
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_data', required=True)
parser.add_argument('-c','--cluster_key', required=True)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

# Parse args
path_data = args['path_data']
path_out = args['path_out']
cluster_key = args['cluster_key']

# Read rna adata
adata = mu.read(path_data)
del adata.mod['atac']
obs = adata.obs.copy()
adata = adata.mod['rna'].copy()
adata.obs = obs

# Clean up
adata.obs[cluster_key] = adata.obs[cluster_key].str.replace("/", "_").copy()
adata.X = adata.layers["counts"].copy()
adata.raw = adata

# Write out
adata.write(path_out)