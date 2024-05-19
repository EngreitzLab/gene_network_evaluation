import pandas as pd
import numpy as np
import mudata as mu

import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_data', required=True)
parser.add_argument('-ed','--path_eReg_direct', required=True)
parser.add_argument('-ee','--path_eReg_extended', required=False, default=None)
parser.add_argument('-r', '--path_r2g', required=True)
parser.add_argument('-g', '--path_grn', required=True)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

# Parse args
path_data = args['path_data']
path_eReg_direct = args['path_eReg_direct']
path_eReg_extended = args['path_eReg_extended']
path_r2g = args['path_r2g']
path_grn = args['path_grn']
path_out = args['path_out']

# Read in mudata
data = mu.read_h5mu(path_data)

# Read in eReg direct
egrn = pd.read_csv(path_eReg_direct, sep="\t")

# Read in eReg extended and concat if necessary
if path_eReg_extended is not None:
    eReg_extended = pd.read_csv(path_eReg_extended, sep="\t")
    egrn = pd.concat([egrn, eReg_extended], axis=0, ignore_index=True)

# Grab r2g
r2g = egrn[["Region", "Gene", "importance_R2G"]]
r2g.columns = ["cre", "gene", "score"]
r2g["pval"] = np.nan
r2g["cluster"] = "global"
r2g.to_csv(path_r2g, index=False)
data.uns["r2g"] = r2g

# Grab grn
grn = egrn[["TF", "Gene", "importance_TF2G"]]
grn.columns = ["tf", "gene", "score"]
grn["pval"] = np.nan
grn["cluster"] = "global"
grn.to_csv(path_grn, index=False)
data.uns["grn"] = grn

# Save data
data.write_h5mu(path_out)

