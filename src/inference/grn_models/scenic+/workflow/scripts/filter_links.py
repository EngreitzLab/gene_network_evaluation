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
parser.add_argument('-c','--path_cistromes', required=True)
parser.add_argument('-g', '--path_grn', required=True)
parser.add_argument('-r', '--path_r2g', required=True)
parser.add_argument('-t', '--path_tri', required=True)
parser.add_argument('-m', '--path_mdata', required=True)

args = vars(parser.parse_args())
path_input = args['path_input']
path_scenicplus_obj = args['path_scenicplus_obj']
path_cistromes = args['path_cistromes']
path_grn = args['path_grn']
path_r2g = args['path_r2g']
path_tri = args['path_tri']
path_mdata = args['path_mdata']

# Read mdata to add objects to
data = mu.read(path_input)

# Load SCENICPLUS object
scplus_obj = dill.load(open(path_scenicplus_obj, 'rb'))

# unfiltered cistromes (TF-region links based on motif enrichment)
cistromes = scplus_obj.uns["Cistromes"]["Unfiltered"]
dill.dump(cistromes, open(path_cistromes, 'wb'))
cistrome_dict = {}
for cistrome in cistromes.keys():
    cistrome_dict[cistrome] = cistromes[cistrome].df[["Chromosome", "Start", "End"]].astype(str).apply(lambda x: "-".join(x), axis=1).values
data.uns["cistromes"] = cistrome_dict

# initial TF-gene grn
grn = scplus_obj.uns["TF2G_adj"][["TF", "target", "importance_x_rho"]]
grn.columns = ["source", "target", "weight"]
grn["pval"] = 1
grn.to_csv(path_grn, sep="\t", index=False)
data.uns["grn"] = grn

# initial region-gene links
r2g = scplus_obj.uns["region_to_gene"] \
    .rename(columns={"importance_x_rho": "weight"})
r2g["pval"] = 1
r2g.drop(columns=["Distance"], inplace=True)
r2g = r2g[["target", "region", "weight", "pval", "importance", "rho"]]
r2g.to_csv(path_r2g, sep="\t", index=False)
data.uns["r2g"] = r2g

# tri -- SCENIC+ eRegulon with TF-gene links pruned by region-gene links and cistromes
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
data.uns["tri"] = tri

# Repeat tri with after standard filtering
apply_std_filtering_to_eRegulons(scplus_obj)
tri_filtered = scplus_obj.uns["eRegulon_metadata"]
tri_filtered["weight"] = tri_filtered[["R2G_importance_x_rho", "TF2G_importance_x_rho"]].mean(axis=1)
tri_filtered["pval"] = 1
tri_filtered["region"] = tri_filtered["Region"].str.replace(":", "-")
tri_filtered = tri_filtered.rename(columns={"TF": "source", "Gene": "target"})
tri_filtered = tri_filtered[["source", "target", "region", "weight", "pval", "is_extended", 
                             "R2G_importance", "R2G_rho", "R2G_importance_x_rho",
                             "R2G_importance_x_abs_rho", "TF2G_importance", "TF2G_regulation",
                             "TF2G_rho", "TF2G_importance_x_abs_rho", "TF2G_importance_x_rho"]]
tri_filtered.to_csv(path_tri.replace(".tsv", "_filtered.tsv"), sep="\t", index=False)
data.uns["tri_filtered"] = tri_filtered

# Save data
data.write_h5mu(path_mdata)
