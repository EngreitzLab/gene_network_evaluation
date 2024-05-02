import pandas as pd
import numpy as np
import muon as mu
import os
import argparse
import anndata


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_input', required=True)
parser.add_argument('-m','--path_modules', required=True)
parser.add_argument('-e','--path_MEs', required=True)
parser.add_argument('-p','--path_mdata', required=True)
args = vars(parser.parse_args())

#
path_input = args['path_input']
path_modules = args['path_modules']
path_MEs = args['path_MEs']
path_mdata = args['path_mdata']

# Load data
mdata = mu.read(path_input)

# Load modules
modules = pd.read_csv(path_modules, sep='\t')

# Load MEs
MEs = pd.read_csv(path_MEs, sep='\t')

# Get programs
programs = pd.DataFrame(index=modules["module"].unique())

# Make Anndata with cell x program matrix
adata = anndata.AnnData(X=MEs.values, obs=pd.DataFrame(index=MEs.index), var=MEs.columns.to_frame(name="color"))

# Get program x feature matrix
loadings = modules.loc[:, modules.columns.str.contains("kME")]
loading_names = [x.split("_")[1] for x in loadings.columns]
loadings.columns = loading_names
loadings = loadings[programs.index].T

# Add program x feature matrix to anndata
adata.varm["loadings"] = loadings.values
adata.uns["loadings"] = {"names": loadings.columns.to_numpy()}

# Store in mdata
mdata.mod["hdwgcna"] = adata

# Save
mdata.write(path_mdata)
