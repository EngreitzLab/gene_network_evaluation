import numpy as np
import mudata as mu
import scipy.sparse as sp
import loompy as lp
import argparse


# Init args
parser = argparse.ArgumentParser(
    prog="python pre.py",
    description="Prepare loom file for SCENIC run from standard MuData."
)
parser.add_argument('-d','--path_data', required=True, type=str, help="path to input MuData object, this is not modified throughout pipeline")
parser.add_argument('-l','--layer', default="counts", required=False, type=str, help="layer in mdata.mod['rna'].layers to use for regression. If not provided, will use counts.")
parser.add_argument('-o','--path_out', required=True, type=str, help="path to output grn.csv containing TF to gene links")
args = vars(parser.parse_args())

# Parse args
path_data = args['path_data']
layer = args['layer']
path_out = args['path_out']

# Load data
mdata = mu.read_h5mu(path_data)

def write_loom(adata, filename, layer=None, use_raw=False):
    """Write AnnData object as loom that is compatible with SCENIC

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix of shape `n_obs` x `n_vars`.
        Rows correspond to cells and columns to genes.
    filename : str
        Filename to save to
    layer : str
        Layer to save instead of `X`. If `None`, `X` is saved.

    Returns
    -------
    None
    """
    print(f"Saving as loom to {filename}")
    if use_raw:
        adata = adata.raw.to_adata()  #only if adata has RAW saved and thats what you want!!
    row_attrs = dict(zip(adata.var.reset_index().columns, adata.var.reset_index().values.T))
    col_attrs = dict(zip(adata.obs.reset_index().columns, adata.obs.reset_index().values.T))
    row_attrs["Gene"] = np.array(adata.var_names)
    col_attrs["CellID"] = np.array(adata.obs_names)
    if layer is not None:
        X = adata.layers[layer]
    else:
        X = adata.X
    if sp.issparse(X):
        X = X.toarray()
    col_attrs["nGene"] = np.array(np.sum(X.transpose() > 0, axis=0)).flatten()
    col_attrs["nUMI"] = np.array(np.sum(X.transpose(), axis=0)).flatten()
    lp.create(filename, X.transpose(), row_attrs, col_attrs)

# Write loom
write_loom(
    mdata.mod["rna"],
    filename=path_out,
    layer="counts"
)
