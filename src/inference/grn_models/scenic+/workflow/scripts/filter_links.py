import os
import dill
import argparse
import muon as mu
import numpy as np
from scenicplus.preprocessing.filtering import apply_std_filtering_to_eRegulons


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_input', required=True)
parser.add_argument('-s', '--path_scenicplus_obj', required=True)
parser.add_argument('-o','--organism', required=True)
parser.add_argument('-g', '--path_grn', required=True)
parser.add_argument('-r', '--path_r2g', required=True)
parser.add_argument('-t', '--path_tri', required=True)

args = vars(parser.parse_args())
path_input = args['path_input']
path_scenicplus_obj = args['path_scenicplus_obj']
organism = args['organism']
path_grn = args['path_grn']
path_r2g = args['path_r2g']
path_tri = args['path_tri']

# Read mdata to add objects to
data = mu.read(path_input)

# Load SCENICPLUS object
scplus_obj = dill.load(open(path_scenicplus_obj, 'rb'))

# Save pipeline outputs
apply_std_filtering_to_eRegulons(scplus_obj)

# grn
grn = scplus_obj.uns["TF2G_adj"][["TF", "target", "importance_x_rho"]]
grn.columns = ["source", "target", "weight"]
grn["pval"] = 1
grn.to_csv(path_grn, sep="\t", index=False)

# r2g
r2g = scplus_obj.uns["region_to_gene"] \
    .rename(columns={"importance_x_rho": "weight"})
r2g["pval"] = 1
r2g.drop(columns=["Distance"], inplace=True)
r2g = r2g[["target", "region", "weight", "pval", "importance", "rho"]]
r2g.to_csv(path_r2g, sep="\t", index=False)

# tri
tri = scplus_obj.uns["eRegulon_metadata"]
tri["weight"] = tri[["R2G_importance_x_rho", "TF2G_importance_x_rho"]].mean(axis=1)
tri["pval"] = 1
tri["region"] = tri["Region"].str.replace(":", "-")
tri = tri.rename(columns={"TF": "source", "Gene": "target"})
tri = tri[["source", "target", "region", "weight", "pval", "is_extended", 
           "R2G_importance", "R2G_rho", "R2G_importance_x_rho",
           "R2G_importance_x_abs_rho", "TF2G_importance", "TF2G_regulation",
           "TF2G_rho", "TF2G_importance_x_abs_rho", "TF2G_importance_x_rho"]]
tri.to_csv(path_tri, sep="\t", index=False)

# Add to mudata and save
data.uns["grn"] = grn
data.uns["r2g"] = r2g
data.uns["tri"] = tri
