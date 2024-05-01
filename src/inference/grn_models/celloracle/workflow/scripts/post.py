import pandas as pd
import muon as mu
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-d','--path_data', required=True)
parser.add_argument('-r','--path_r2g', required=True)
parser.add_argument('-t','--path_tf2r', required=True)
parser.add_argument('-g','--path_grn', required=True)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

# Parse args
path_input = args['path_data']
path_r2g = args['path_r2g']
path_tf2r = args['path_tf2r']
path_grn = args['path_grn']
path_mdata = args['path_out']

# Read mdata to add objects to
data = mu.read(path_input)

# Region-gene links
r2g = pd.read_csv(path_r2g)
data.uns["r2g"] = r2g

# TF-region links
motifs_df = pd.read_csv(path_tf2r)
data.uns["tf2r"] = motifs_df

# Filter links and concatenate
grn = pd.read_csv(path_grn)
data.uns["grn"] = grn

# TODO: Decide if you want to post process: e.g. extract triplets or filter by grn

# Save data
data.write_h5mu(path_mdata)
